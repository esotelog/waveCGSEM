/* ------------------------------------------------------------------------
 *
 * SPDX-License-Identifier: LGPL-2.1-or-later
 * Copyright (C) 2006 - 2024 by the deal.II authors
 *
 * This file is part of the deal.II library.
 *
 * Part of the source code is dual licensed under Apache-2.0 WITH
 * LLVM-exception OR LGPL-2.1-or-later. Detailed license information
 * governing the source code and code contributions can be found in
 * LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
 *
 * ------------------------------------------------------------------------
 *
 * Author: Wolfgang Bangerth, Texas A&M University, 2006
 */

#ifndef WAVE_CGSEM_H       
#define WAVE_CGSEM_H


#include <filesystem>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <exception>
#include <string>
#include <cmath>

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/symmetric_tensor.h>
#include <deal.II/base/exceptions.h>       // for Assert and ExcMessage
#include <deal.II/base/point.h>  
#include <deal.II/base/array_view.h>


#include <deal.II/lac/vector.h>
#include <deal.II/lac/sparse_matrix.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>


#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_tools.h>

#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/mapping_q1.h>

#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_creator.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/base/utilities.h>

#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_in.h>
 #include <deal.II/grid/grid_tools_geometry.h>
 
#include "functions.h"
//#include "Param_handler.h"


namespace dii = dealii;


template <int dim>
  class waveCGSEM 
  {
  public:
    waveCGSEM(const Parameters & prm);
    void run();
  private:
    void read_mesh();
    void create_mesh();
    void setup_system();
    void assemble_MKmatrices();
    void assemble_boundary_matrix();
    void assemble_point_source();
    void Neumann_plane_wave();
    void timestep_loop();
    void output_results() const;
    void check_conditioning(const dii::SparseMatrix<double> &A, 
                            const std::string &name);
    void write_receiver_header(std::ofstream &receiver_out);
    void write_receiver_data(const double time, std::ofstream &receiver_out);

    dii::Triangulation<dim> mesh;
    dii::FESystem<dim> cg_fe;
    dii::DoFHandler<dim>  dof_handler;
    dii::MappingQ1<dim> mapping;
    
    dii::AffineConstraints<double> constraints;
    const dii::QGaussLobatto<dim> qformula;
    const dii::QGaussLobatto<dim - 1> qfaceformula;
    dii::SparsityPattern      sparsity_pattern;
    dii::SparseMatrix<double> mass_matrix;
    dii::SparseMatrix<double> stiffness_matrix;
    dii::SparseMatrix<double> boundary_matrix;
    dii::Vector<double> solution_u, solution_v;
    dii::Vector<double> old_solution_u, old_solution_v;
    dii::Vector<double> mass_vector;
    dii::Vector<double> inverse_mass_vector;
    dii::Vector<double> source_term;
   
  
    Parameters param;
    Source<dim> source_function;
    typename dii::DoFHandler< dim>::active_cell_iterator source_cell;
    dii::Point<dim> source_pos_ref;
    bool source_initialized = false;
    int timestep_number;
    std::vector<dii::Point<dim>> receivers_pos;
    Mesh_operations<dim> mesh_op;
    const unsigned int bid_min=31, bid_max=32;

    const unsigned int boundary_id_left = 0;
    const unsigned int boundary_id_right = 1;

    std::ofstream out;

  };


#include "template_waveCGSEM.tcc"

#endif