
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
    //
    if (param.source == 4 )
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


    dii::FEFaceValues<dim> fe_facevalues (cg_fe, qfaceformula,
                                           dii::update_values | 
                                           dii::update_quadrature_points |
                                           dii::update_JxW_values);

    const unsigned int dofs_per_cell = cg_fe.n_dofs_per_cell();
    dii::Vector<double> cell_source_vector(dofs_per_cell);
    std::vector<dii::types::global_dof_index> local_dof_indices(dofs_per_cell);
    const dii::Tensor<1,dim> force_vector = source_function.get_force_vector();

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

        }
  }
 
template <int dim>
   void waveCGSEM<dim>::assemble_boundary_matrix()
   {
    const double rho=param.rho;
    const double vp= param.vp;
    const double vs= param.vs;

    
    const bool is_boundary_source = (param.source ==4);
    const bool point_source_abc =(!is_boundary_source);
      
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
               bool is_periodic = (bid == boundary_id_left || bid == boundary_id_right);
               bool boundary_source_abc= (is_boundary_source && !is_periodic);

               if (point_source_abc || boundary_source_abc)
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
   if (param.source < 4 )//point sourc
       assemble_point_source ();
   else if (param.source == 4)//Neumann boundary source
        Neumann_plane_wave ();
   else
      throw std::runtime_error("Invalid source type");

                     
    // Open output file for receivers' data (all y-levels in one file)
   const std::string base = param.outputfile;
   std::ofstream receiver_out((base + ".dat").c_str());
   write_receiver_header(receiver_out);

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
    
    while (true)
      {
        const double time = static_cast<double>(timestep_number - 1) * time_step;
        if (time > total_time)
           break;
        
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


          //***Solving for u**
        solution_u = old_solution_u;
        solution_u.add(time_step, solution_v);
        constraints.distribute(solution_u);

         //calculating energ
        dealii::Vector<double> tmp(solution_v); 
        tmp.scale(mass_vector); 
        double energy = 0.5*(solution_v * tmp) + 
                         0.5*stiffness_matrix.matrix_norm_square(solution_u);

        std::cout << "Energy: " << energy 
                   <<std::endl<<std::endl;
        out << "Energy: " << energy 
              <<std::endl<<std::endl;

        if (energy>1e6) //stop the simulation if the total energy is too high (1e6 is a arbitrary value )
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
        write_receiver_data(time, receiver_out);
          //if (timestep_number >= 3)
           //  break;
        ++timestep_number;
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
      std::filesystem::create_directories(output_dir);
      const std::string filename = output_dir +
                                   param.outputfile + "-" +
                                   dii::Utilities::int_to_string (timestep_number, 3) +
                                   ".vtk";
      std::ofstream output (filename.c_str());
      data_out.write_vtk (output);
    }
  

  template<int dim>
  void waveCGSEM<dim>::write_receiver_header(std::ofstream &receiver_out)
  {
    const unsigned int nx = param.receivers_x.size();
    const unsigned int ny = param.receivers_y.size();
    if (receivers_pos.empty() || nx == 0 || ny == 0 ||
        receivers_pos.size() != nx * ny)
      return;

    const char *comp_name[] = {"ux", "uy", "uz"};
    const unsigned int nr = receivers_pos.size();

    receiver_out << "#";
    for (unsigned int k = 0; k < nr; ++k)
      {
        const auto &p = receivers_pos[k];
        receiver_out << " x" << k << "=" << p[0] << ", y" << k << "=" << p[1];
        if (k + 1 < nr)
          receiver_out << ",";
      }
    receiver_out << std::endl;

    receiver_out << "# time";
    for (unsigned int k = 0; k < nr; ++k)
      for (unsigned int c = 0; c < dim; ++c)
        receiver_out << "  " << comp_name[c] << "_r" << k;
    receiver_out << std::endl;
  }

template<int dim>
  void waveCGSEM<dim>::write_receiver_data(const double time,
                                         std::ofstream &receiver_out)
  {
    const unsigned int nx = param.receivers_x.size();
    const unsigned int ny = param.receivers_y.size();
    if (receivers_pos.empty() || nx == 0 || ny == 0 ||
        receivers_pos.size() != nx * ny)
      return;

    // Same point order as read_receiver_coordinates (y outer, x inner).

    dii::Vector<double> value(dim);

    receiver_out << std::scientific << std::setprecision(8);
    receiver_out << time;

    for (unsigned int k = 0; k < receivers_pos.size(); ++k)
      {
        dii::VectorTools::point_value(
          dof_handler, solution_u, receivers_pos[k], value);
        for (unsigned int c = 0; c < dim; ++c)
          receiver_out << " " << value(c);
      }
    receiver_out << std::endl;
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

