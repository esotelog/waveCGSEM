#ifndef FUNCTIONS_H       
#define FUNCTIONS_H

#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <exception>
#include <string>
#include <vector>
#include <array>
#include <cmath>
#include <algorithm>
#include <limits>
#include <stdexcept>

#include <deal.II/base/utilities.h>
#include <deal.II/base/function.h>
#include <deal.II/base/conditional_ostream.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/exceptions.h> 
#include <deal.II/base/exception_macros.h>        // for Assert and ExcMessage
#include <deal.II/base/point.h>  
#include <deal.II/base/numbers.h>
#include <deal.II/base/tensor.h>
#include <deal.II/base/quadrature_lib.h>

#include <deal.II/lac/full_matrix.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_interface_values.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_refinement.h>
#include <deal.II/grid/grid_out.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/meshworker/dof_info.h>
#include <deal.II/meshworker/integration_info.h>
#include <deal.II/meshworker/mesh_loop.h>
#include <deal.II/meshworker/local_integrator.h>
#include <deal.II/meshworker/assembler.h>
#include <deal.II/meshworker/simple.h>
#include <deal.II/meshworker/copy_data.h>
#include <deal.II/meshworker/scratch_data.h>

#include "Param_handler.h"
namespace dii = dealii;
namespace diiM = dii::MeshWorker;


template <int dim>
  class Source: public dii::Function <dim>
  {
   private:
    Parameters param;
    const double fo;
    const double source_rad;
    const int stf_type;
    const double M11, M12, M22, M21;
    const double loc_x, loc_y, loc_z;
    
  public:
    Source(const Parameters & param_);
    double source_duration() const;
    double source_time (const double & t) const;
    double apply_source_time() const;
    dii::SymmetricTensor<2, dim> get_MT() const;
    dii::Tensor<1, dim> get_force_vector() const;
    bool is_isotropic() const;
    static bool is_isotropic(const Parameters & prm);
    static dii::Point<dim> get_source_position(const Parameters & prm);
    dii::Point<dim> source_pos;
   
  };

template <int dim>
class Mesh_operations
{
private:
  std::vector<dii::Point<dim>> get_boundary_points() const;
  double get_velocity() const;
  const Parameters param;
  const double vs, vp;
  const int wavesamples, sourcesampling;
  const double fo, source_rad;
  const std::vector<dii::Point<dim>> boundary_points;
  double hbg, outer_radius;

public:
  Mesh_operations(const Parameters &param_);
  static double convert_coord(const Parameters &param,
                             const int &coord,
                             const double &coord_value);
  void mesh_create(dii::Triangulation<dim> &mesh) const;
  void label_extreme_boundary_faces(dii::DoFHandler<dim> &dof_handler,
                                    const unsigned int boundary_id, const unsigned int bid_min,
                                    const unsigned int bid_max, const unsigned int axis, // 0 for x, 1 for y
                                    const double fraction_nf, double &min, double &Lmax,
                                    std::vector<std::pair<typename dii::DoFHandler<dim>::cell_iterator,
                                    unsigned int>> &faces) const;
  void source_refine(dii::Triangulation<dim> &mesh);
  void print_mesh(dii::Triangulation<dim> &mesh) const;
  std::vector<dii::Point<dim>> receiver_coordinates(const double &h,
                                                    const double &yfrac,
                                                    const double &top_dg_y,
                                                    const double &bottom_dg_y,
                                                    const unsigned int &ncells,
                                                    const std::vector<double> &receivers_x) const;

  const double Lx, Ly;
};

template <int dim>
  class InitialValuesU : public dii::Function<dim>
  {
    
  public:
    InitialValuesU();
    virtual void vector_value(const dii::Point<dim> & /*p*/,
                              dii::Vector<double> & value) const override;
  };

template <int dim>
  class InitialValuesV : public dii::Function<dim>
  {
  public:
    InitialValuesV();
    virtual void vector_value(const dii::Point<dim> & /*p*/,
                              dii::Vector<double> & value) const override;
  };

template <int dim>
class BoundaryValuesBase : public dii::Function<dim>
{
protected:
    const Parameters param;
    Source<dim> source_function;
    double taper_function(const double &x) const;
public:
    BoundaryValuesBase(const Parameters & param_);
};

 template <int dim>
  class BoundaryValuesU : public BoundaryValuesBase<dim>
  {
    public:
       BoundaryValuesU(const Parameters & param_);
       virtual void vector_value(const dii::Point<dim> & p,
                              dii::Vector<double> & value) const override;
  };


  template <int dim>
  class BoundaryValuesV : public BoundaryValuesBase<dim>
  {
    private:
      double derivative_source_time(const double & t) const;
    public:
       BoundaryValuesV(const Parameters & param_);
       virtual void vector_value(const dii::Point<dim> & p,
                              dii::Vector<double> & value) const override;
  };

template <int dim>
dii::SymmetricTensor<4, dim> get_stiffness_tensor(const double & vp,
                                                  const double & vs,
                                                  const double & rho);


template <int dim>
dii::SymmetricTensor<2, dim> get_strain(const dii::FEValues<dim> & fe_values,
                                   const unsigned int & shape_func,
                                   const unsigned int  & q_point);
template <int dim>
dii::SymmetricTensor<2, dim> get_strain_comp(const dii::FEInterfaceValues<dim> & face_values,
                                             const bool & here_there,
                                             const unsigned int & shape_index,
                                             const unsigned int  & q_point);
template <int dim>
dii::SymmetricTensor<2, dim> get_strain_comp(const dii::FEFaceValuesBase<dim> & face_values,
                                             const unsigned int & shape_index,
                                             const unsigned int  & q_point);

template <int dim>
dii::SymmetricTensor<2, dim> get_fracture_stiffness(const double & Zn,
                                                    const double & Zt);
template<int dim>
void cg_dg_cells(dii::DoFHandler<dim> &dof_handler, 
                 double & fcoord,
                 const int & rows_above,
                 const int & cg_id,
                 const int & dg1_id, 
                 const int & dg2_id,
                 double & top_dg_y,
                 double & bottom_dg_y,
                 double & h);    
template<int dim> 
std::vector<dii::Point<dim>> read_receiver_coordinates(const std::vector<double> & receivers_x,
                                                        const std::vector<double> & receivers_y);


double calc_time_interval(const double & space,
                         const double & vel,
                         const unsigned int & p_degree);

void rounding_digits(double & value, 
                     const int & digits, 
                     bool roundown_flag=true);

void echo_parameters(const Parameters &params);


#include "template_functions.tcc"

#endif