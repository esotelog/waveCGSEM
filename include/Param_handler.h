#ifndef PARAM
#define PARAM

#include <deal.II/base/parameter_handler.h>
#include <deal.II/base/enable_observer_pointer.h>
#include <fstream>
#include <iostream>

namespace dii = dealii;

struct Parameters
{
  public:
    Parameters();
    double vp;
    double vs;
    double rho;

    //source parameters
    int stf_type;
    double fo;
    double source_rad;
    int source;
    double M11;
    double M12;
    double M22;
    double M21;
    double phi;
    double taper;

    double loc_x;
    double loc_y;
    double loc_z;// in 2D this is not used

    // time discretization
 
    double total_time; // (tfactor>0; total time = source duration)
    double dt_factor;
    double theta;
    //mesh paramenters
    int nsamples; /*number of samples per min wavelength*/
    double Lx;
    double Ly;
    int sourcesampling; /* diameter of source will be sampled 6 times*/
    bool read_mesh;

    //fem
    int dim;
    int p_degree;
  
    
   //output 
   //receiver x coordinates
   std::vector<double> receivers_x;
   std::vector<double> receivers_y;

   bool vtk;
   int vtk_step;

   //file name (vtk and trace files)
   std::string outputfile;
};


class ParameterReader: public  dii::Subscriptor
{
public:

  ParameterReader(const std::string & parameter_file_);
  Parameters get_parameters();


private:

  void declare_parameters();
  void read_parameters();
  void pass_parameters ();

  const std::string parameter_file;

  //default constructed objects. Note: default constructors should be defined for these class objects
  dii::ParameterHandler  prm;
  Parameters param;

};

#endif