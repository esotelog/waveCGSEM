template <int dim>
  Source<dim>::Source(const Parameters  & param_ ):
  dii::Function<dim>(dim),
  param(param_),
  fo(param_.fo),
  source_rad(param_.source_rad),
  stf_type(param_.stf_type),
  M11(param_.M11), M12(param_.M12), M22(param_.M22), M21(param_.M21),
  loc_x(param_.loc_x), loc_y(param_.loc_y), loc_z(param_.loc_z),
  source_pos(get_source_position(param))
  { 

  }

  template <int dim>
  dii::Point<dim> Source<dim>::get_source_position( const Parameters & prm) 
  {
    dii::Point<dim> p;

    if constexpr (dim == 2)
      {
        p[0] = prm.loc_x;
        p[1] = prm.loc_y;
      }
    else if constexpr (dim == 3)
      {
        p[0] = prm.loc_x;
        p[1] = prm.loc_y;
        p[2] = prm.loc_z;
      }
    else
      {
        static_assert(dim == 2 || dim == 3, "Unsupported dimension");
      }

    return p;

  }

template <int dim>
  dii::SymmetricTensor<2, dim> Source<dim>::get_MT() const
  {
     dii::SymmetricTensor<2, dim> tmp;
     tmp.clear();  
         // Diagonal terms
          // Diagonals
    tmp[0][0] = M11;
    tmp[1][1] = M22;

    // Off-diagonals
    if (param.source == 1)
      {
        tmp[0][1] = M12;
        tmp[1][0] = M12;  // MUST BE EXPLICIT
      }
    else if (param.source == 3)
      {
        tmp[0][1] = M12;
        tmp[1][0] = -M12;  // MUST BE EXPLICIT
      }
    else
      {
        throw std::runtime_error("Unsupported source type");
      }
    return tmp;

  }

  template <int dim>
  dii::Tensor<1, dim> Source<dim>::get_force_vector() const
  {
    dii::Tensor<1, dim> tmp;
    const double pi = dii::numbers::PI;
    const double phi_rad = param.phi * pi / 180.0;
    tmp[0] = std::cos(phi_rad);
    tmp[1] = std::sin(phi_rad);
    return tmp;
  }

template <int dim>
  bool Source<dim>::is_isotropic() const
  {
    bool is_isotropic = (param.source == 1 && (param.M11 == param.M22 && param.M12 == 0.0));
    bool is_Patboundary=(param.source == 4 && param.phi == -90.0);

    return (is_isotropic || is_Patboundary);
  }

template <int dim>
  bool Source<dim>::is_isotropic(const Parameters & prm)
  {
    bool is_isotropic = (prm.source == 1 && (prm.M11 == prm.M22 && prm.M12 == 0.0));
    bool is_Patboundary=(prm.source == 4 && prm.phi == -90.0);

    return (is_isotropic || is_Patboundary);
  }

template <int dim>
    double Source<dim>::source_duration() const
    {
        return 4.0/fo;
    }




template <int dim>
    double Source<dim>::source_time(const double & t) const
    { 
       const double to = 2.0/fo;
       const double pi = dii::numbers::PI;
       const double f2 = pi*pi*fo*fo;
       const double time2 = (t-to)*(t-to);
       
       const double stf = (1.0-2.0*f2*time2)*std::exp(-f2*time2);

       const double stf_value = (t>source_duration() ) ? 0.0 : stf;
       return stf_value;
    }

  template <int dim>
    double Source<dim>::apply_source_time() const
    {
        return source_time(this->get_time());
    }
    

template <int dim>
   Mesh_operations<dim>::Mesh_operations(const Parameters & param_):
    param(param_),
    vs (param_.vs), 
    vp(param_.vp),
    wavesamples(param_.nsamples), 
    sourcesampling(param_.sourcesampling), 
    fo(param_.fo), 
    source_rad(param_.source_rad),
    boundary_points(get_boundary_points()),
    hbg(0.0),
    outer_radius(0.0),
    Lx(param_.Lx),
    Ly(param_.Ly)
    { 
    }
    

  template <int dim>
   std::vector<dii::Point<dim>> Mesh_operations<dim>::get_boundary_points() const 
   { 
    std::vector<dii::Point<dim>> bpoints(2);
    const double Lx = param.Lx;
    const double Ly = param.Ly;

     if constexpr ( dim==2)
      { 
        bpoints[0]=dii::Point<dim>(0.0, 0.0);
        bpoints[1]=dii::Point<dim>(Lx, -Ly);
      }
     
    else if constexpr ( dim==3)
      { 
        bpoints[0]=dii::Point<dim>(0.0, 0.0, 0.0);
        bpoints[1]=dii::Point<dim>(Lx, Ly, -Ly);
      
      }

    else
      {
         static_assert(dim == 2 || dim == 3, "Unsupported dimension.");
      }
    return bpoints;
   }

  template <int dim>
   double Mesh_operations<dim>::convert_coord( const Parameters & param,
                                               const int & coord, 
                                               const double &coord_value)
   {

    double tmp_Lx=coord_value*param.Lx;
    double tmp_Ly=coord_value*param.Ly;

    if (coord == 0)
        return tmp_Lx;

    if constexpr (dim == 2)
    {
        if (coord == 1 || coord == -1)
            return -tmp_Ly;
    }
    else if constexpr (dim == 3)
    {
        if (coord == 1)
            return tmp_Ly;
        if (coord == 2 || coord == -1)
            return -tmp_Ly;
    }
    else
    {
        static_assert(dim == 2 || dim == 3, "Unsupported dimension");
    }

    throw std::runtime_error("Unsupported coordinate value");
  
   }

   template <int dim>
   double Mesh_operations<dim>::get_velocity() const
   {
    bool isotropic= Source<dim>::is_isotropic(param);
    if (isotropic)
      return vp;
    else return vs;
   }
   

template <int dim>
   void Mesh_operations<dim>::mesh_create(dii::Triangulation<dim> & mesh) const 
   {  
     double ratio, Lmax, Lmin, eps=1e-10;
     enum direction {x,y};
     direction axis;

      if (Lx>Ly+ eps || std::abs(Lx-Ly) < eps)
        {
          Lmax = Lx;
          Lmin = Ly;
          axis = x;

        }
     else
        {
          Lmax = Ly;
          Lmin = Lx;
          axis = y;
        }

     ratio = Lmax/Lmin;
     
     // create a new mesh and refine it globally
     const double c = get_velocity();
     const double lambda_s = c/fo;
     const double h = lambda_s/wavesamples;

     if (Lmin/h > 1)
     {
       const unsigned int refine_min = static_cast<unsigned int> ( std::ceil(std::log2(Lmin/h)) );
       const unsigned int rep_min = 1u << refine_min;  //2^refine_min
       const unsigned int rep_max =
            std::max(1u, static_cast<unsigned int>(std::ceil(rep_min * (ratio))));

        std::vector<unsigned int> repetitions_dim(dim);
        repetitions_dim[0] = (axis == x) ? rep_max : rep_min;
        repetitions_dim[1] = (axis == x) ? rep_min : rep_max;


        mesh.clear();
        dii::GridGenerator::subdivided_hyper_rectangle(mesh,
                                                      repetitions_dim,
                                                      boundary_points[0], 
                                                      boundary_points[1],
                                                      true) ;
                                         

       }  
     // If ratio > 1, refine the mesh anisotropically
     /*
     if (ratio > 1)
       {
         const unsigned int n_aniso_refine=ratio-1;

         auto cut_axis = (axis == x) ? 
              dii::RefinementCase<dim>::cut_axis(x) : dii::RefinementCase<dim>::cut_axis(y);

         for (unsigned int i = 0; i < n_aniso_refine; ++i)
         {
          for (auto &cell : mesh.active_cell_iterators())
               cell->set_refine_flag(cut_axis);
  
          mesh.execute_coarsening_and_refinement();
         }

       }
    */
     
   }

  template <int dim>
   void Mesh_operations<dim>::label_extreme_boundary_faces(dii::DoFHandler<dim> &dof_handler,
                                                           const unsigned int boundary_id,
                                                           const unsigned int bid_min,
                                                           const unsigned int bid_max,
                                                           const unsigned int axis, // 0 for x, 1 for y
                                                           const double fraction_nf,
                                                           double &Lmin, double &Lmax,
                                                           std::vector<std::pair<typename dii::DoFHandler<dim>::cell_iterator,
                                                           unsigned int>> &faces) const
   {


    if (axis >= dim)
       throw std::runtime_error("axis must be in [0, " + std::to_string(dim - 1) +
                   "] for dim=" + std::to_string(dim) + ".");

    if (fraction_nf >= 0.5)
        throw std::runtime_error(
           "fraction_nf must be < 0.5 (otherwise min/max boundary bands can overlap).");

    faces.clear();
 
  // --- Step 1: collect faces with given boundary_id ---
   for (auto &cell : dof_handler.active_cell_iterators())
     {
     for (unsigned int f = 0; f < dii::GeometryInfo<dim>::faces_per_cell; ++f)
       {
         if (cell->face(f)->at_boundary() && cell->face(f)->boundary_id() == boundary_id)
         {
             faces.emplace_back(cell, f);
         }
       }
    }
 
   const unsigned int nfaces = faces.size();
   if (nfaces == 0)
        throw std::runtime_error("No faces found with boundary_id: " + 
                                 std::to_string(boundary_id) + ".");
 

       // --- Step 2: sort faces by m-th coordinate ---
   std::sort(faces.begin(), faces.end(),[axis](const auto &a, const auto &b)
          {
             const double xa = a.first->face(a.second)->center()[axis];
             const double xb = b.first->face(b.second)->center()[axis];
              return xa < xb;
           });
   
   // --- Step 3 : assign faces with bid_min→ 31 ---
   const double Lref=param.Lx*fraction_nf;
   Lmin = 0.0;

   //for (unsigned int i = 0; i < nextreme_out; ++i)
   for ( auto & it: faces)
       {
         auto &cell = it.first;
         const unsigned int f = it.second;
         auto face = cell->face(f);
   
         face->set_boundary_id(bid_min);
   
         // accumulate length along m-direction
         Lmin += cell->extent_in_direction(axis);
         if (Lmin >= Lref)
           break;
       }
   
       // --- Step 4: assign faces with bid_max→ 32 ---

   Lmax = 0.0;

   //for (unsigned int i = 0; i < nextreme_out; ++i)
   for (auto rit =faces.rbegin(); rit != faces.rend(); ++rit)  
       {
        auto &cell =rit->first;
        const unsigned int f = rit->second;
        auto face = cell->face(f);
   
        face->set_boundary_id(bid_max);
   
        Lmax += cell->extent_in_direction(axis);
        if (Lmax >= Lref)
          break;
       }
   }
 

  template <int dim>
   void Mesh_operations<dim>::source_refine(dii::Triangulation<dim> & mesh) 
   {
    const dii::Point<dim> source_pos = Source<dim>::get_source_position(param);
    const double target_h = 2.0 * source_rad / sourcesampling;

    // Background cell size from the already-refined uniform mesh
    const double h_bg = mesh.begin_active()->diameter();

    if (h_bg <= target_h * 1.05)
    {
        hbg = h_bg;
        std::cout << "Source refinement: h_bg (" << h_bg
                  << ") <= target_h (" << target_h << "), no refinement needed." << std::endl;
        return;
    }

    // Doubling layers needed: target_h -> 2*target_h -> ... -> h_bg
    const int nlayers_auto = static_cast<int>(std::ceil(std::log2(h_bg / target_h)));

    const double safety_factor = 1.2; 
    const double inner_radius = safety_factor * source_rad;

    // Transition zone: layer k has cell size target_h*2^k, spans ~1.5x its size
    const double layer_thickness_factor = 1.5;
    const double transition = layer_thickness_factor * target_h
                              * (std::pow(2.0, nlayers_auto) - 1.0);
    outer_radius = inner_radius + transition;

    std::cout << "Source refinement: target_h=" << target_h 
              << ", h_bg=" << h_bg
              << ", nlayers=" << nlayers_auto
              << ", inner_r=" << inner_radius 
              << ", outer_r=" << outer_radius << std::endl;

    auto get_target_h_at_distance = 
    [target_h, h_bg, inner_radius, this]
    (double d) -> double 
    {
        if (d <= inner_radius)
            return target_h;

        if (d >= outer_radius)
            return h_bg;

        const double t = (d - inner_radius) / (outer_radius - inner_radius);
        return target_h * std::exp(t * std::log(h_bg / target_h));
    };

    while (true)
    {
        bool flag_refine = false;
        for (auto &cell : mesh.active_cell_iterators())
        {
            double h = cell->diameter();
            double d = cell->center().distance(source_pos);
            double target_h_local = get_target_h_at_distance(d);

            if (h > target_h_local * 1.05)
            {
                cell->set_refine_flag();
                flag_refine = true;
            }
        }

        if (flag_refine)
            mesh.execute_coarsening_and_refinement();
        else
            break;
    }
    
    for (const auto &cell : mesh.active_cell_iterators())
    {
        double d = cell->center().distance(source_pos);
        if (d >= 1.1 * outer_radius)
        {
            hbg = cell->extent_in_direction(1);
            break;
        }
    }
    
    if (hbg == 0.0)
         throw std::runtime_error(
         "No uniform mesh region found outside refinement zone. "
         "Ensure mesh domain is large enough.");

   }

template <int dim>
void Mesh_operations<dim>::print_mesh(dii::Triangulation<dim> &mesh) const
{
  const std::string filename = "./vtk/" + param.outputfile + "_mesh.vtk";
  std::ofstream out(filename);
  dii::GridOut grid_out;
  grid_out.write_vtk(mesh, out);
}

template <int dim>
  std::vector<dii::Point<dim>> Mesh_operations<dim>::receiver_coordinates(const double & h,
                                                                          const double & yfrac, 
                                                                          const double & top_dg_y,
                                                                          const double & bottom_dg_y,
                                                                          const unsigned int & ncells,
                                                                          const std::vector<double> & receivers_x
                                                                          ) const
  { 
    std::vector<dii::Point<dim>> receivers_pos;
    if (receivers_x.empty())
      return receivers_pos;

    // Convert x coordinates from normalized to physical and round to 1 decimal
    const int x_coord = 0;
    std::vector<double> xcoords(receivers_x.size());
    for (unsigned int i = 0; i < receivers_x.size(); ++i)
    {
      xcoords[i] = convert_coord( param, x_coord, receivers_x[i]);
      xcoords[i] = std::round(xcoords[i] * 10.0) / 10.0;
    }

    // y-coordinates for each receiver line, rounded to 1 decimal
    double yup   = top_dg_y + (ncells + 0.5) * h;
    double ydown = bottom_dg_y - (ncells + 0.5) * h;
    double ydg   = yfrac + 0.5 * h;
    yup   = std::round(yup   * 10.0) / 10.0;
    ydown = std::round(ydown * 10.0) / 10.0;
    ydg   = std::round(ydg   * 10.0) / 10.0;

    auto point_at = [](double x, double y) -> dii::Point<dim>
    {
      if constexpr (dim == 2)
        return dii::Point<dim>(x, y);
      else
        return dii::Point<dim>(x, y, 0.0);
    };

    // CG receivers at yup
    for (const double x : xcoords)
      receivers_pos.push_back(point_at(x, yup));

    // CG receivers at ydown
    for (const double x : xcoords)
      receivers_pos.push_back(point_at(x, ydown));

    // DG receivers at ydg
    for (const double x : xcoords)
      receivers_pos.push_back(point_at(x, ydg));

    return receivers_pos;
  }

template <int dim>
  InitialValuesU<dim>::InitialValuesU():
  dii::Function<dim>(dim) // considering FE_Nothing components
  {}

template <int dim>
  void InitialValuesU<dim>::vector_value(const dii::Point<dim> & /*p*/,
                                        dii::Vector<double> &values) const 
  {
    AssertDimension(values.size(),dim);
    values = 0.0;

  }

template <int dim>
  InitialValuesV<dim>::InitialValuesV():
  dii::Function<dim>(dim)
  {}

template <int dim>
  void InitialValuesV<dim>::vector_value(const dii::Point<dim> & /*p*/,
                                        dii::Vector<double> &values) const 
  {
    AssertDimension(values.size(), dim);
    values = 0.0;

  }

  template <int dim>
  BoundaryValuesBase<dim>::BoundaryValuesBase(const Parameters & param_):
  dii::Function<dim>(dim), // considering FE_Nothing components
  param(param_),
  source_function(param_)
  {}

 template <int dim>
  double BoundaryValuesBase<dim>::taper_function(const double & x)const
  {
    const double L = param.Lx;
    const double Lsupport = param.taper * L;
    const double pi=dii::numbers::PI;
  

    if (x < Lsupport)
      return 0.5*(1.0-std::cos(pi*x/Lsupport));
    else if (x > L-Lsupport)
      return 0.5*(1.0-std::cos(pi*(L-x)/Lsupport));
    else
      return 1.0;
  }

template <int dim>
 BoundaryValuesU<dim>::BoundaryValuesU(const Parameters & param_):
 BoundaryValuesBase<dim>(param_)
  {}

template <int dim>
  void BoundaryValuesU<dim>::vector_value(const dii::Point<dim> & p,
                                        dii::Vector<double> &values) const
  {
    AssertDimension(values.size(), dim);
    double taper = BoundaryValuesBase<dim>::taper_function(p[0]);

    dii::Tensor<1, dim> u_vector =
          BoundaryValuesBase<dim>::source_function.get_force_vector();
    
    double stf = BoundaryValuesBase<dim>::source_function.source_time(this->get_time());

    values[0] = taper * u_vector[0] * stf;
    values[1] = taper* u_vector[1] * stf;

  }

template <int dim>
  BoundaryValuesV<dim>::BoundaryValuesV(const Parameters & param_):
  BoundaryValuesBase<dim>(param_)
  {}

template <int dim>
  double BoundaryValuesV<dim>::derivative_source_time(const double & t) const
  {
    const double fo = BoundaryValuesBase<dim>::param.fo;
    const double to = 2.0/fo;
    const double pi = dii::numbers::PI;
    const double a = pi*fo*(t-to);
    const double a2=a*a;
    
    double stf= 2.0*pi*fo*a*(2.0*a2-3.0)*std::exp(-a2);
    const double sd = BoundaryValuesBase<dim>::source_function.source_duration();
    double stf_value = (t>sd ) ? 0.0 : stf;
    return stf_value;
  }


template <int dim>
  void BoundaryValuesV<dim>::vector_value(const dii::Point<dim> & p,
                                        dii::Vector<double> &values) const
  {
    AssertDimension(values.size(), dim);
    double taper = BoundaryValuesBase<dim>::taper_function(p[0]);
    dii::Tensor<1, dim> u_vector=BoundaryValuesBase<dim>::source_function.get_force_vector();
    double stf=derivative_source_time(this->get_time());

    values[0] = taper * u_vector[0] * stf;
    values[1] = taper* u_vector[1] * stf;

  }


template <int dim>
dii::SymmetricTensor<4, dim> get_stiffness_tensor(const double & vp,
                                                  const double & vs,
                                                  const double & rho)
 {
   const double mu =vs*vs*rho;
   const double lambda =rho*vp*vp - 2*mu;

   dii::SymmetricTensor<4, dim> tmp;
      for (unsigned int i = 0; i < dim; ++i)
        for (unsigned int j = 0; j < dim; ++j)
          for (unsigned int k = 0; k < dim; ++k)
            for (unsigned int l = 0; l < dim; ++l)
              tmp[i][j][k][l] = (((i == k) && (j == l) ? mu : 0.0) +
                                 ((i == l) && (j == k) ? mu : 0.0) +
                                 ((i == j) && (k == l) ? lambda : 0.0));
      return tmp;
 }

template <int dim>
dii::SymmetricTensor<2, dim> get_strain(const dii::FEValues<dim> & fe_values,
                                        const unsigned int & shape_func,
                                        const unsigned int  & q_point)
  {
    dii::SymmetricTensor<2, dim> tmp;
    
    // Get the actual component this shape function belongs to
    const unsigned int comp = fe_values.get_fe().system_to_component_index(shape_func).first;
    
    // Get the full gradient (all spatial derivatives)
    const dii::Tensor<1, dim> grad_phi = fe_values.shape_grad(shape_func, q_point);
    
    // Diagonal strain: ε_ii = ∂u_i/∂x_i
    tmp[comp][comp] = grad_phi[comp];
    
    // Off-diagonal strains: ε_ij = 0.5 * (∂u_i/∂x_j + ∂u_j/∂x_i)
    for (unsigned int i = 0; i < dim; ++i)
    {
      if (i != comp)
        tmp[comp][i] = tmp[i][comp] = 0.5 * grad_phi[i];
    }
           
    return tmp;
  }

  template <int dim>
  dii::SymmetricTensor<2, dim> get_strain_comp(const dii::FEInterfaceValues<dim> &face_values,
                                               const bool & here_there,
                                               const unsigned int & shape_index,
                                               const unsigned int  & q_point,
                                               const unsigned int comp)
  {
    dii::SymmetricTensor<2,dim> eps;

    // Gradient of the comp component of the shape function
    const dii::Tensor<1,dim> grad_phi =
        face_values.shape_grad(here_there, shape_index, q_point, comp);

   //The following maps any comp component  to 0...dim-1
   // This is relevant for the dg_fe_system, where DG components are numbered after dim-1:, dim,...2*dim-1
   const unsigned int mcomp = comp % dim;

    // Fill eps tensor
   eps[mcomp][mcomp] = grad_phi[mcomp];

   for (unsigned int i = 0; i < dim; ++i)
    {
     if (i != mcomp)
        eps[mcomp][i] = eps[i][mcomp] = 0.5 * grad_phi[i];
   }

  return eps;
}


template <int dim>
dii::SymmetricTensor<2, dim> get_strain_comp(const dii::FEFaceValuesBase<dim> & face_values,
                                             const unsigned int & shape_index,
                                             const unsigned int  & q_point)
  {
    dii::SymmetricTensor<2,dim> eps;
    const unsigned int comp = face_values.get_fe().system_to_component_index(shape_index).first;
    const dii::Tensor<1,dim> grad_phi =face_values.shape_grad(shape_index, q_point);

   //The following maps any comp component  to 0...dim-1
   // This is relevant for the dg_fe_system, where DG components are numbered after dim-1:, dim,...2*dim-1
   const unsigned int mcomp = comp % dim;

   // Fill eps tensor
   eps[mcomp][mcomp] = grad_phi[mcomp];

   for (unsigned int i = 0; i < dim; ++i)
    {
     if (i != mcomp)
        eps[mcomp][i] = eps[i][mcomp] = 0.5 * grad_phi[i];
   }


    return eps;
  }

template <int dim>
  dii::SymmetricTensor<2, dim> get_fracture_stiffness(const double & Zn,
                                                      const double & Zt)
{
  const double eps = 1e-30;
  const double inv_Zn = 1.0 / (std::abs(Zn) + eps);
  const double inv_Zt = 1.0 / (std::abs(Zt) + eps);

  dii::SymmetricTensor<2, dim> tmp = 0.0;
  if constexpr (dim == 2)
  {
    tmp[0][0] = inv_Zt;
    tmp[1][1] = inv_Zn;
  }
  else
  {
    tmp[0][0] = inv_Zt;
    tmp[1][1] = inv_Zt;
    tmp[2][2] = inv_Zn;
  }
  return tmp;
}

template<int dim>
void cg_dg_cells(dii::DoFHandler<dim> &dof_handler, 
                 double & fcoord,
                 const int & rows_above,
                 const int & cg_id,
                 const int & dg1_id, 
                 const int & dg2_id,
                 double & top_dg_y,
                 double & bottom_dg_y,
                 double & h)
{
    // Step 1: Find cell closest to fcoord to get h in the uniform region
    double min_dist = std::numeric_limits<double>::max();
    typename dii::DoFHandler<dim>::active_cell_iterator fracture_cell;

    for (const auto &cell : dof_handler.active_cell_iterators())
    {
      const double dist = std::abs(cell->center()[1] - fcoord);
      if (dist < min_dist)
      {
        min_dist = dist;
        fracture_cell = cell;
      }
    }

    h = fracture_cell->extent_in_direction(1);

    // Step 2: Find the horizontal interior face closest to fcoord
    double best_dist = std::numeric_limits<double>::max();
    double face_y = fcoord;

    for (const auto &cell : dof_handler.active_cell_iterators())
      for (unsigned int f = 0; f < dii::GeometryInfo<dim>::faces_per_cell; ++f)
      {
        if (cell->face(f)->at_boundary())
          continue;
        // Skip vertical faces: only keep faces whose normal is predominantly in y
        const dii::Tensor<1, dim> normal = cell->face(f)->center() - cell->center();
        if (std::abs(normal[1]) < std::abs(normal[0]))
          continue;
        const double fy = cell->face(f)->center()[1];
        const double dist = std::abs(fy - fcoord);
        if (dist < best_dist)
        {
          best_dist = dist;
          face_y = fy;
        }
      }

    fcoord = face_y;

    // Step 3: DG zone symmetric around fcoord
    top_dg_y    = fcoord + rows_above * h;
    bottom_dg_y = fcoord - rows_above * h;
    const double tol = 1e-10;

    // Step 4: Classify cells
    int n_cg = 0, n_dg1 = 0, n_dg2 = 0;

    for (auto &cell : dof_handler.active_cell_iterators())
    {
      const double cy = cell->center()[1];

      if (cy >= bottom_dg_y - tol && cy <= top_dg_y + tol)
      {
        cell->set_active_fe_index(1);
        if (cy >= fcoord)
          { cell->set_material_id(dg1_id); ++n_dg1; }
        else
          { cell->set_material_id(dg2_id); ++n_dg2; }
      }
      else
      {
        cell->set_active_fe_index(0);
        cell->set_material_id(cg_id);
        ++n_cg;
      }
    }

    std::cout << "cg_dg_cells: fcoord = " << fcoord
              << ", h = " << h
              << ", top_dg_y = " << top_dg_y
              << ", bottom_dg_y = " << bottom_dg_y
              << ", CG: " << n_cg 
              << ", DG1: " << n_dg1 
              << ", DG2: " << n_dg2 << std::endl;
}

template<int dim> 
 std::vector<dii::Point<dim>> read_receiver_coordinates (
                                                        const std::vector<double> & receivers_x,
                                                        const std::vector<double> & receivers_y
                                                        )
{
  std::vector<dii::Point<dim>> receivers_pos;

  for (auto & y: receivers_y)
    for (auto & x : receivers_x) 
    { 
    if constexpr (dim == 2)
       receivers_pos.push_back(dii::Point<dim>(x, y));
    else if constexpr (dim == 3)
    receivers_pos.push_back(dii::Point<dim>(x, y, 0.0));
    }
  return receivers_pos;
  }

