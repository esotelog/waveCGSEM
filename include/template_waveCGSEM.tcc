
//constructor 
template<int dim>
  waveCGSEM<dim>::waveCGSEM(const Parameters &prm):
    cg_fe(dii::FE_Q<dim>(prm.p_degree)^ dim),
    dof_handler(mesh),
    qformula(prm.p_degree + 1),
    qfaceformula(prm.p_degree + 1),
    param(prm),
    source_function(prm),
    timestep_number(1),
    mesh_op(prm)
    
   { 
    out.open ("waveCGSEM_log.txt");
    std::cout << " Source position: " << source_function.source_pos << std::endl;
    out << " Source position: " << source_function.source_pos << std::endl;
   }

template<int dim>
  void waveCGSEM<dim>::read_mesh()
  {
    const std::string filename = "./mesh_cg.vtk";
    std::ifstream input(filename);
    if (!input.is_open())
      throw std::runtime_error("Cannot open mesh file: " + filename);

    mesh.clear();
    dii::GridIn<dim> grid_in;
    grid_in.attach_triangulation(mesh);
    grid_in.read_vtk(input);
    receivers_pos = read_receiver_coordinates<dim>(param.receivers_x, param.receivers_y);

    std::cout << "Mesh read from " << filename << std::endl;
    std::cout << "Number of active cells: " << mesh.n_active_cells() << std::endl;

    std::cout << "Receivers' coordinates:" << std::endl;
    for (unsigned int i = 0; i < receivers_pos.size(); ++i)
      std::cout << "  [" << i << "] " << receivers_pos[i] << std::endl;

    out << "Mesh read from " << filename << std::endl;
    out << "Number of active cells: " << mesh.n_active_cells() << std::endl;
    out << "Receivers' coordinates:" << std::endl;
    for (unsigned int i = 0; i < receivers_pos.size(); ++i)
      out << "  [" << i << "] " << receivers_pos[i] << std::endl;
  }

template<int dim>
  void waveCGSEM<dim>::create_mesh()
  {
    //Mesh_operations<dim> mesh_op(param);
    mesh_op.mesh_create(mesh);
    const double h_bg = mesh.begin_active()->diameter();

    //mesh_op.source_refine(mesh);
    mesh_op.print_mesh(mesh);

    // Collect periodic faces

    // reading receiver coordinates
    receivers_pos = read_receiver_coordinates<dim>(param.receivers_x, param.receivers_y);

    std::cout << "Mesh Lx: " << mesh_op.Lx <<std::endl;
    std::cout << "Mesh Ly: " << mesh_op.Ly <<std::endl; 
    std::cout << "Background cell size: " << h_bg << std::endl;

    std::cout << "Receivers' coordinates:" << std::endl;
    for (unsigned int i = 0; i < receivers_pos.size(); ++i)
      std::cout << "  [" << i << "] " << receivers_pos[i] << std::endl;

    std::cout << "Number of active cells: " << mesh.n_active_cells()
                << std::endl;
   

    // writing to file
    out << "Mesh Lx: " << mesh_op.Lx <<std::endl;
    out << "Mesh Ly: " << mesh_op.Ly <<std::endl; 
    out << "Background cell size: " << h_bg << std::endl;
            
    out << "Receivers' coordinates:" << std::endl;
    for (unsigned int i = 0; i < receivers_pos.size(); ++i)
        out << "  [" << i << "] " << receivers_pos[i] << std::endl;
               
    out << "Number of active cells: " << mesh.n_active_cells()
                            << std::endl;

  }

template<int dim>
  void waveCGSEM<dim>::setup_system()
  {    

    dof_handler.distribute_dofs(cg_fe);
  
    std::cout << "Number of degrees of freedom: " << dof_handler.n_dofs()
                << std::endl
                << std::endl;

  
  // peridoic boundary condition

  //periodic boundary parameters

   const unsigned int direction = 0; // 0 for x, 1 for y
   std::vector<dii::GridTools::PeriodicFacePair<
    typename dii::DoFHandler<dim>::cell_iterator>> periodic_faces;

   dii:: GridTools::collect_periodic_faces(dof_handler,
                                          boundary_id_left,
                                          boundary_id_right,
                                          direction,   // 0 for x, 1 for y
                                          periodic_faces);
    
    //Apply periodic BC
    constraints.clear();
    dii::DoFTools::make_periodicity_constraints<dim,dim>(periodic_faces,
                                                         constraints);
    constraints.close();          
    
    //Sparsity pattern
    dii::DynamicSparsityPattern dsp(dof_handler.n_dofs(), dof_handler.n_dofs());
    dii::DoFTools::make_sparsity_pattern(dof_handler, dsp, constraints,false);
    sparsity_pattern.copy_from(dsp);
 
    //initializing global matrices
    stiffness_matrix.reinit(sparsity_pattern);
    boundary_matrix.reinit(sparsity_pattern);
    // initializing  global vectors
    mass_vector.reinit(dof_handler.n_dofs());
    inverse_mass_vector.reinit(dof_handler.n_dofs());
    solution_u.reinit(dof_handler.n_dofs());
    solution_v.reinit(dof_handler.n_dofs());
    old_solution_u.reinit(dof_handler.n_dofs());
    old_solution_v.reinit(dof_handler.n_dofs());
    source_term.reinit(dof_handler.n_dofs());
  

    //source_cells = source_function.find_source_cells(dof_handler);

    out << "Number of degrees of freedom: " << dof_handler.n_dofs()
                << std::endl
                << std::endl;

  }

template<int dim>
  void waveCGSEM<dim>::assemble_MKmatrices()
  {
    std::cout << " assembly begins " << std::endl;
    out << " assembly begins " << std::endl;

    //Assembly of stiffness matrix
    /**************************** */
    
    const dii::SymmetricTensor<4, dim> stiffness_tensor = 
                                       get_stiffness_tensor<dim>(param.vp,param.vs,param.rho);
    dii::FEValues<dim> fe_values(cg_fe,
                                 qformula,
                                 dii::update_values | dii::update_gradients |
                                 dii::update_quadrature_points | dii::update_JxW_values);

    const unsigned int dofs_per_cell = cg_fe.n_dofs_per_cell();
    dii::FullMatrix<double> cell_stiffness_matrix(dofs_per_cell, dofs_per_cell);
    dii::Vector<double>     cell_mass_vector(dofs_per_cell);
    std::vector<dii::types::global_dof_index> local_dof_indices(dofs_per_cell);

    //assembly per cell
    for (const auto &cell : dof_handler.active_cell_iterators())    
        {
          cell_stiffness_matrix = 0;
          cell_mass_vector = 0;
          fe_values.reinit(cell);

          for (const unsigned int q_point: fe_values.quadrature_point_indices())
           {  // Precompute strains
             std::vector<dii::SymmetricTensor<2, dim>> eps_phi(dofs_per_cell);
             for (unsigned int i=0; i<dofs_per_cell; ++i)
               eps_phi[i] = get_strain(fe_values, i, q_point);

             for (const unsigned int i :fe_values.dof_indices())
              {
               double phi_i = fe_values.shape_value(i,q_point);
               
               cell_mass_vector(i) +=  phi_i * phi_i *fe_values.JxW(q_point);

                for (const unsigned int j :fe_values.dof_indices())           
                  cell_stiffness_matrix(i, j) += (eps_phi[i] *  stiffness_tensor *  eps_phi[j]) *                    
                                                  fe_values.JxW(q_point); 
                  // mass matrix                
              }
            }//end of loop by quadrature points 
          //Trasfer local cell matrix to global laplace_matrix
           cell->get_dof_indices(local_dof_indices);
           constraints.distribute_local_to_global(cell_stiffness_matrix, cell_mass_vector,
                                                local_dof_indices, stiffness_matrix, mass_vector);
            /*
            for (const unsigned int i: fe_values.dof_indices())  
              { 
                mass_vector(local_dof_indices[i]) += cell_mass_vector(i);

                for (const unsigned int j: fe_values.dof_indices())           
                    stiffness_matrix.add(local_dof_indices[i],
                                         local_dof_indices[j],
                                         cell_stiffness_matrix(i,j));
                   
               }
            */
          
          }//enf of loop by cell dof   
    // inverse mass vector
    const double mass_eps =
    1.0e-20 * std::max(mass_vector.linfty_norm(), 1.0);
  
    for (unsigned int i = 0; i < dof_handler.n_dofs(); ++i)
       inverse_mass_vector(i) =
      (mass_vector(i) > mass_eps) ? (1.0 / mass_vector(i)) : 0.0;
  
    constraints.distribute(inverse_mass_vector);
    
    std::cout << "assembly completed " <<std::endl <<std::endl;
    out << "assembly completed" << std::endl;
  }

  template<int dim>
  void waveCGSEM<dim>::assemble_point_source()
  {
    //Mapping source to reference cell
    if (!source_initialized) //only once
    {
      dii::Point<dim> p = source_function.source_pos;
      
      auto result =
      dii::GridTools::find_active_cell_around_point(mapping,
                                               dof_handler,
                                               p,
                                               {},
                                               1e-12);
     
      source_initialized = true;
      source_cell = result.first;
      source_pos_ref=result.second;

      std::cout << "p_ref: " << source_pos_ref << std::endl;
      out << "p_ref: " << source_pos_ref << std::endl;
      // Print all vertices in real (physical) coordinates
      for (unsigned int v = 0; v < dii::GeometryInfo<dim>::vertices_per_cell; ++v)
       {
         const dii::Point<dim> &vertex = source_cell->vertex(v);
         std::cout << "vertex " << v << ": " << vertex << std::endl;
         out << "vertex " << v << ": " << vertex << std::endl;
       }
      }

    // 1) Build a "quadrature" with exactly the reference point
    std::vector<dii::Point<dim>> quad_points(1);
    quad_points[0] = source_pos_ref;
    dii::Quadrature<dim> point_quadrature(quad_points);

    // 2) FEValues to evaluate shape function gradients at that point
    dii::FEValues<dim> fe_point_values(mapping,
                                      cg_fe,
                                      point_quadrature,
                                      dii::update_gradients);
    fe_point_values.reinit(source_cell);

    const unsigned int dofs_per_cell = cg_fe.n_dofs_per_cell();
    std::vector<dii::types::global_dof_index> local_dof_indices(dofs_per_cell);
    source_cell->get_dof_indices(local_dof_indices);



    for (unsigned int i = 0; i < dofs_per_cell; ++i)
    {
      const unsigned int comp_i = cg_fe.system_to_component_index(i).first;  
      double force_i = 0.0;
      if (param.source == 1 || param.source == 3)
      {
        dii::SymmetricTensor<2,dim> M = source_function.get_MT();  
        dii::Tensor<1,dim> grad_phi_i = fe_point_values.shape_grad_component(i, 0, comp_i);
        
        for (unsigned int j = 0; j < dim; ++j)
          force_i += M[comp_i][j] * grad_phi_i[j]; 
      }

      else if (param.source == 2)
      {
             // component of dof
        const dii::Tensor<1,dim> force_vector = source_function.get_force_vector();
        const double phi = cg_fe.shape_value_component(i, source_pos_ref, comp_i);
        force_i = force_vector[comp_i] * phi;
      }
      else
        throw std::runtime_error("Invalid source type");

      source_term(local_dof_indices[i]) +=force_i;
     /*
      std::cout <<"phi: " << phi << std::endl;
      std::cout <<"force_vector: " << force_vector[comp_i] << std::endl;
      std::cout <<"source_value: " << source_value << std::endl;
      std::cout <<"source_term: " << source_term(local_dof_indices[i]) << std::endl;
      std::cout <<"--------------------------------" << std::endl;

      out <<"phi: " << phi << std::endl;
      out <<"force_vector: " << force_vector[comp_i] << std::endl;
      out <<"source_value: " << source_value << std::endl;
      out <<"source_term: " << source_term(local_dof_indices[i]) << std::endl;
      out <<"--------------------------------" << std::endl;
     */
    }

  }


template<int dim>
  void waveCGSEM<dim>::Neumann_plane_wave()
  {
    //labeling the extreme boundary faces with bid_min and bid_max to taper the source
    //boundary_cellsf: vector of pairs of cell iterators and face indices at boundary_id
    std::vector<std::pair<typename dii::DoFHandler<dim>::cell_iterator,
    unsigned int>> boundary_cellsf;
    const unsigned int boundary_id = 3, axis = 0;
    double Lmin = 0.0, Lmax = 0.0;
    const double fraction_nf = param.taper;
   
    mesh_op.label_extreme_boundary_faces(dof_handler, 
                                         boundary_id, bid_min, bid_max,
                                         axis, fraction_nf, Lmin, Lmax, boundary_cellsf);

    //source taper function
    /*
    auto source_taper=[&](const double &x, const unsigned int &bid)
      {
        const double L=param.Lx;
        double s = 0.0;
       
        if (bid == bid_min)  s = x/Lmin;
        else if (bid == bid_max) s = (L-x)/Lmax;
        else return 1.0;  
        
        s = std::min(1.0, std::max(0.0, s));   // clamp s
        return 10*std::pow(s, 3) - 15*std::pow(s, 4) + 6*std::pow(s, 5);
      
     };
     */

    dii::FEFaceValues<dim> fe_facevalues (cg_fe, qfaceformula,
                                           dii::update_values | 
                                           dii::update_quadrature_points |
                                           dii::update_JxW_values);

    const unsigned int dofs_per_cell = cg_fe.n_dofs_per_cell();
    dii::Vector<double> cell_source_vector(dofs_per_cell);
    std::vector<dii::types::global_dof_index> local_dof_indices(dofs_per_cell);
    const dii::Tensor<1,dim> force_vector = source_function.get_force_vector();
        // impedence scaling for source at boundary
      /* 
    auto impedance =[&](const dii::Tensor<1,dim> & f)
      {
       const double rho=param.rho;
       const double vp= param.vp;
       const double vs= param.vs;

       const double norm = std::sqrt(f[0]*f[0] + f[1]*f[1]);
       if (norm < 1e-12) throw std::runtime_error("Zero force vector");
   
       return rho * std::sqrt( vp*vp*(f[0]*f[0])/ (norm*norm) + 
                     vs*vs*(f[1]*f[1])/ (norm*norm) );
      };

      const dii::Tensor<1,dim> traction=impedance(force_vector)*force_vector;
      */
    for (const auto &it : boundary_cellsf)
 
        {
          const auto &cell = it.first;
          const unsigned int f = it.second;
          //const unsigned int bid = cell->face(f)->boundary_id();

          cell_source_vector = 0;
          fe_facevalues.reinit(cell, f);       
          for (const unsigned int q_point : fe_facevalues.quadrature_point_indices())
            {
              //const double x = fe_facevalues.quadrature_point(q_point)[axis];
              //const double taper = source_taper(x, bid);
             for (const unsigned int i: fe_facevalues.dof_indices())
              {
                const unsigned int comp_i = cg_fe.system_to_component_index(i).first;
                
                double phi_i = fe_facevalues.shape_value(i,q_point);
                cell_source_vector(i) += phi_i* 
                                         force_vector[comp_i]*
                                         fe_facevalues.JxW(q_point);
              }
            }
            cell->get_dof_indices(local_dof_indices);

            constraints.distribute_local_to_global(cell_source_vector,
                                                  local_dof_indices, source_term);
            /*
            for (const unsigned int i: fe_facevalues.dof_indices())
              source_term(local_dof_indices[i]) += cell_source_vector(i);
            */
        }
  }

template<int dim>
  void waveCGSEM<dim>::Dirichlet_plane_wave(dii::Vector<double> &solvec,
                                            const double &time,
                                            dii::Function<dim> &bfunction)
  {
    const unsigned int boundary_id = 3;
    bfunction.set_time(time);
    std::map<dii::types::global_dof_index, double> boundary_values;

    dii::VectorTools::interpolate_boundary_values(dof_handler,
                                                  boundary_id,
                                                  bfunction,
                                                  boundary_values);

     for (const auto &it : boundary_values)
          solvec[it.first] = it.second;
  }
 
template <int dim>
   void waveCGSEM<dim>::assemble_boundary_matrix()
   {
    const double rho=param.rho;
    const double vp= param.vp;
    const double vs= param.vs;

    //const bool is_boundary_source =(param.source == 4 || param.source == 5);
   //const bool is_point_source =(! is_boundary_source);
      
    dii::FEFaceValues<dim> fe_facevalues (cg_fe, qfaceformula,
                                          dii::update_values | 
                                          dii::update_quadrature_points |
                                          dii::update_normal_vectors |
                                          dii::update_JxW_values);

    const unsigned int dofs_per_cell = cg_fe.n_dofs_per_cell();
    dii::FullMatrix<double> cell_boundary_matrix(dofs_per_cell, dofs_per_cell);

    std::vector<dii::types::global_dof_index> local_dof_indices(dofs_per_cell);

    for (const auto &cell : dof_handler.active_cell_iterators())
        for (unsigned int f = 0; f < dii::GeometryInfo<dim>::faces_per_cell; ++f)
            if (cell->at_boundary(f))            
              {                 
               const unsigned int bid = cell->face(f)->boundary_id();
               //bool is_periodic = (bid == boundary_id_left || bid == boundary_id_right);
                bool is_top_boundary = (bid == 3 || bid == bid_max || bid == bid_min);
               //bool allow_boundary_source= (is_boundary_source && !is_top_boundary);
                
                //if (is_point_source || allow_boundary_source)
               if (bid==2 || is_top_boundary)
                {
                cell_boundary_matrix = 0;
                fe_facevalues.reinit(cell, f);
                const std::vector<dii::Tensor<1, dim>> &normals = fe_facevalues.get_normal_vectors();

                std::vector<unsigned int> component(dofs_per_cell);
                for (unsigned int i=0;i<dofs_per_cell;++i)
                     component[i] = cg_fe.system_to_component_index(i).first;

                for (const unsigned int q_point : fe_facevalues.quadrature_point_indices())
                 {
                   const dii::Tensor<1, dim> &normal = normals[q_point];

                   for (const unsigned int i: fe_facevalues.dof_indices()) 
                    {
                      const unsigned int comp_i = component[i];
                      double ni = normal[comp_i];
                      double phi_i = fe_facevalues.shape_value(i,q_point);

                      for (const unsigned int j: fe_facevalues.dof_indices())
                          {
                            const unsigned int comp_j = component[j];
                            double nj = normal[comp_j];
                            double phi_j = fe_facevalues.shape_value(j,q_point);
                            // creating a Kronecker delta
                            double delta_ij = (comp_i == comp_j) ;
                           

                            double geometry_term = rho*(vp - vs)*ni *nj  + rho*vs*delta_ij;

                            cell_boundary_matrix(i, j) += phi_i * phi_j * 
                                                          geometry_term *
                                                          fe_facevalues.JxW(q_point) ;
                          }
                    }
                  }
                cell->get_dof_indices(local_dof_indices);
                constraints.distribute_local_to_global(cell_boundary_matrix,
                                                      local_dof_indices, boundary_matrix);
                /*
                for (const unsigned int i: fe_facevalues.dof_indices())         
                  for (const unsigned int j: fe_facevalues.dof_indices())
                       boundary_matrix.add(local_dof_indices[i], local_dof_indices[j], cell_boundary_matrix(i, j));
                */
              }
            }
    }
  

  template <int dim>
   void waveCGSEM<dim>::timestep_loop()
  {

    const double rho=param.rho;
    const double dt_factor=param.dt_factor;

    //dii::AffineConstraints<double> constraints;
    //constraints.close();
  
     //intial values 
    dii::VectorTools::project (dof_handler,
                               constraints,
                               qformula,
                               InitialValuesU<dim>(),
                               old_solution_u);

    dii::VectorTools::project (dof_handler,
                              constraints,
                              qformula,
                              InitialValuesV<dim>(),                       
                              old_solution_v);

    // constraining old solution vectors
    constraints.distribute(old_solution_u);
    constraints.distribute(old_solution_v);

   dii::Vector<double> tmp (solution_u.size());
   dii::Vector<double> forcing_terms (solution_u.size());

   source_term = 0.0;
   if (param.source != 4 && param.source != 5)//point sourc
       assemble_point_source ();
   else if (param.source == 4)//Neumann boundary source
        Neumann_plane_wave ();
   else if (param.source == 5)//Dirichlet boundary source
      {
        bvalues_vfunction = std::make_unique<BoundaryValuesV<dim>>(param);
        bvalues_ufunction = std::make_unique<BoundaryValuesU<dim>>(param);
      };
                     
    // Open output files for receivers' data (one per y-level)
   const std::string base = param.outputfile + "-";
   std::ofstream out_yup((base + "rcg_yup.dat").c_str());
   std::ofstream out_ydg((base + "rcg_ydg.dat").c_str());
   std::ofstream out_ydown((base + "rcg_ydown.dat").c_str());
   write_receiver_header(out_yup, out_ydg, out_ydown);

    // time calculations
    // const double min_cellsize= extreme_cellsize(mesh);
    const double min_cellsize = dii::GridTools::minimal_cell_diameter(mesh);
    std::cout << "hmin = " << min_cellsize <<std::endl;
    double time_step = dt_factor * calc_time_interval(min_cellsize, param.vp, param.p_degree);
    const double total_time=param.total_time;

    std::cout <<"time_step:  " <<time_step << std::endl;
    std::cout << "   Total_time: " << total_time << std::endl<< std::endl;

    out << "hmin = " << min_cellsize <<std::endl;
    out <<"time_step:  " <<time_step << std::endl;
    out << "   Total_time: " << total_time << std::endl<< std::endl;
    
    for (double time=0.0;
         time<=total_time;
         time+=time_step, ++timestep_number)
         {
   
            std::cout <<"timestep_number:  " <<timestep_number << std::endl;
            std::cout <<"current time:  " <<time << std::endl;

            out <<"timestep_number:  " <<timestep_number << std::endl;
            out <<"current time:  " <<time << std::endl;

         solution_v = old_solution_v;

         stiffness_matrix.vmult (tmp, old_solution_u);
         tmp *= -time_step/rho;
         tmp.scale(inverse_mass_vector);
         solution_v += tmp;

         //Absorbing  Boundary condition for v
         boundary_matrix.vmult (tmp, old_solution_v);
         tmp *= -time_step/rho;
         tmp.scale(inverse_mass_vector);
         solution_v += tmp;

          //Include the source
          if (time <= source_function.source_duration() && param.source != 5)
          {
            source_function.set_time(time);
            double source_value= source_function.apply_source_time();

            forcing_terms.equ(source_value * time_step/rho, source_term);; //forcing terms = dt* S(time)
            
            forcing_terms.scale(inverse_mass_vector);
          }
          else
            forcing_terms =0.0;

          solution_v += forcing_terms;
          constraints.distribute(solution_v);

          
         //**Adding Dirichlet boundary condition for v**
         if (param.source ==5)
            Dirichlet_plane_wave(solution_v, time, *bvalues_vfunction);


          //***Solving for u**
          solution_u = old_solution_u;
          solution_u.add(time_step, solution_v);
          constraints.distribute(solution_u);
          //**Adding Dirichlet boundary condition for u**
          if (param.source ==5)
             Dirichlet_plane_wave(solution_u, time, *bvalues_ufunction);

         //calculating energ
         dealii::Vector<double> tmp(solution_v); 
         tmp.scale(mass_vector); 
         double energy = 0.5*(solution_v * tmp) + 
                         0.5*stiffness_matrix.matrix_norm_square(solution_u);

         std::cout << "Energy: " << energy 
                   <<std::endl<<std::endl;
         out << "Energy: " << energy 
              <<std::endl<<std::endl;

         if (
          energy>1e6) //stop the simulation if the total energy is too high (1e6 is a arbitrary value )
              {
               std::cout << "solution does not converge: stopping the simulation"
               << std::endl;
               out << "solution does not converge: stopping the simulation"
               << std::endl;
               break;
              }

         if (timestep_number % param.vtk_step == 0 && param.vtk)
             output_results ();
      
          old_solution_u = solution_u;
          old_solution_v = solution_v;
          write_receiver_data(time, out_yup, out_ydg, out_ydown);
          //if (timestep_number >= 3)
           //  break;
          
         }
  }



  template <int dim>
    void  waveCGSEM<dim>::output_results()  const
    {
      dii::DataOut<dim> data_out;
      data_out.set_flags (dii::DataOutBase::VtkFlags(param.p_degree)); // Since fe.degree is 2
      data_out.attach_dof_handler (dof_handler);
      data_out.add_data_vector (solution_u, "U");
      //data_out.add_data_vector (solution_v, "V");

        // --- Norm of U ---
      // Use FE's component_to_system_index() as per deal.II manual (fe_system.h):
      // "you should not rely on this numbering... Rather use system_to_component_index()
      // and component_to_system_index()"
      const unsigned int n_dofs = solution_u.size();
      dii::Vector<double> U_norm(n_dofs);
      const unsigned int n_dofs_per_component = cg_fe.n_dofs_per_cell() / dim;
      std::vector<dii::types::global_dof_index> local_dof_indices(cg_fe.n_dofs_per_cell());

      for (const auto &cell : dof_handler.active_cell_iterators())
        if (!cell->is_artificial())
          {
            cell->get_dof_indices(local_dof_indices);
            for (unsigned int idx = 0; idx < n_dofs_per_component; ++idx)
              {
                double sum2 = 0.0;
                for (unsigned int comp = 0; comp < dim; ++comp)
                  {
                    const unsigned int local_i =
                      cg_fe.component_to_system_index(comp, idx);
                    const auto global_i = local_dof_indices[local_i];
                    sum2 += solution_u(global_i) * solution_u(global_i);
                  }
                const double mag = std::sqrt(sum2);
                for (unsigned int comp = 0; comp < dim; ++comp)
                  {
                    const unsigned int local_i =
                      cg_fe.component_to_system_index(comp, idx);
                    const auto global_i = local_dof_indices[local_i];
                    U_norm(global_i) = mag;
                  }
              }
          }
      data_out.add_data_vector(U_norm, "U_norm");

      // Subdivide patches to show smooth Q2 solution (default 0 gives jagged cell-boundary look)
      data_out.build_patches (param.p_degree);
      
      const std::string output_dir = "./vtk/";
      const std::string filename = output_dir +
                                   param.outputfile + "-" +
                                   dii::Utilities::int_to_string (timestep_number, 3) +
                                   ".vtk";
      std::ofstream output (filename.c_str());
      data_out.write_vtk (output);
    }
  

  template<int dim>
  void waveCGSEM<dim>::write_receiver_header(std::ofstream &out_yup,
                                            std::ofstream &out_ydg,
                                            std::ofstream &out_ydown)
  {
    const unsigned int nr = param.receivers_x.size();
    if (receivers_pos.empty() || nr == 0 || receivers_pos.size() != 3 * nr)
      return;

    const char *comp_name[] = {"ux", "uy", "uz"};
    auto write_header = [&](std::ofstream &f, unsigned int offset, const char *label)
    {
      f << "# " << label << "  y = " << receivers_pos[offset][1] << std::endl;
      f << "# x_coords:";
      for (unsigned int i = 0; i < nr; ++i)
        f << " " << receivers_pos[offset + i][0];
      f << std::endl;
      f << "# time";
      for (unsigned int i = 0; i < nr; ++i)
        for (unsigned int c = 0; c < dim; ++c)
          f << "  " << comp_name[c] << "_r" << i;
      f << std::endl;
    };

    write_header(out_yup,      0,      "CG yup");
    write_header(out_ydg, nr,     "DG ydg");
    write_header(out_ydown,   2 * nr, "CG ydown");
  }

template<int dim>
  void waveCGSEM<dim>::write_receiver_data(const double time,
                                         std::ofstream &out_yup,
                                         std::ofstream &out_ydg,
                                         std::ofstream &out_ydown)
  {
    const unsigned int nr = param.receivers_x.size();
    if (receivers_pos.empty() || nr == 0)
      return;

    // receivers_pos layout:
    //   [0    .. nr-1]    : receivers at y-level 0
    //   [nr   .. 2*nr-1]  : receivers at y-level 1
    //   [2*nr .. 3*nr-1]  : receivers at y-level 2
    // Vector has dim components (ux, uy for 2D; ux, uy, uz for 3D). Same for all receivers.

    dii::Vector<double> value(dim);

    // y-level 0 (yup): time ux_r0 uy_r0 ux_r1 uy_r1 ...
    out_yup << time;
    for (unsigned int i = 0; i < nr; ++i)
    {
      dii::VectorTools::point_value(dof_handler, solution_u, receivers_pos[i], value);
      for (unsigned int c = 0; c < dim; ++c)
        out_yup << " " << value(c);
    }
    out_yup << std::endl;

    // y-level 1 (ydg): time ux_r0 uy_r0 ux_r1 uy_r1 ...
    out_ydg << time;
    for (unsigned int i = nr; i < 2 * nr; ++i)
    {
      dii::VectorTools::point_value(dof_handler, solution_u, receivers_pos[i], value);
      for (unsigned int c = 0; c < dim; ++c)
        out_ydg << " " << value(c);
    }
    out_ydg << std::endl;

    // y-level 2 (ydown): time ux_r0 uy_r0 ux_r1 uy_r1 ...
    out_ydown << time;
    for (unsigned int i = 2 * nr; i < 3 * nr; ++i)
    {
      dii::VectorTools::point_value(dof_handler, solution_u, receivers_pos[i], value);
      for (unsigned int c = 0; c < dim; ++c)
        out_ydown << " " << value(c);
    }
    out_ydown << std::endl;
  }

  template<int dim>
  void waveCGSEM<dim>::check_conditioning(const dii::SparseMatrix<double> &A,
                                        const std::string &name)
  {
    std::cout << "Estimating conditioning for " << name << " (power iterations)..." << std::endl;
    
    dii::Vector<double> v(dof_handler.n_dofs());
    for (unsigned int i = 0; i < v.size(); ++i)
      v[i] = std::sin(0.37 * (i + 1));
    v /= v.l2_norm();

    // Power iteration for lambda_max (5 iterations)
    double lambda_max = 0.0;
    for (int iter = 0; iter < 5; ++iter) {
      dii::Vector<double> Av(v.size());
      A.vmult(Av, v);
      lambda_max = Av.l2_norm();
      v = Av;
      v /= lambda_max;
    }

    // Inverse power for lambda_min (5 iterations with diagonal preconditioner as cheap inverse)
    v = 1.0;  // reset to constant
    v /= v.l2_norm();
    
    dii::Vector<double> diag(dof_handler.n_dofs());
    diag = 1.0;  // simplified: use identity for cheap diagonal inverse
    
    double lambda_min_inv = 0.0;
    for (int iter = 0; iter < 5; ++iter) {
      dii::Vector<double> tmp(v.size());
      A.vmult(tmp, v);
      dii::Vector<double> Dinv_Av(v.size());
      for (unsigned int i = 0; i < v.size(); ++i)
        Dinv_Av[i] = diag[i] * tmp[i];
      lambda_min_inv = Dinv_Av.l2_norm();
      v = Dinv_Av;
      v /= lambda_min_inv;
    }
    double lambda_min = 1.0 / (lambda_min_inv + 1e-16);

    const double cond_est = lambda_max / (lambda_min + 1e-16);

    std::cout << "  λ_max~ " << lambda_max << ", λ_min~ " << lambda_min 
              << ", κ~ " << cond_est << std::endl;
    out << "Conditioning (" << name << "): κ~ " << cond_est << std::endl;
  }
template<int dim>
  void waveCGSEM<dim>::run()
  { 
    if (param.read_mesh)
      read_mesh();
    else
      create_mesh();
    setup_system();
    assemble_MKmatrices();
    assemble_boundary_matrix();
    timestep_loop();
  }

