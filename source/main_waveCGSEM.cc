#include "waveCGSEM.h"
int main(int argc, char *argv[])
{  try
    {   
      if (argc < 2)
        {
        std::cerr << "Usage: " << argv[0] << " <parameter-file.prm>\n";
        return 1;
        }
        std::string parameter_file = argv[1];
      
        ParameterReader read_param (parameter_file);
        Parameters param_init=read_param.get_parameters();
        int dim = param_init.dim;

        switch(dim)
        {
          case 2:
          {
            waveCGSEM<2> wave_sim(param_init);
            wave_sim.run();
            break;
          }
          case 3:
          {
            waveCGSEM<3> wave_sim(param_init);
            wave_sim.run();
            break;
          }
        }

    }   
   catch (std::exception &exc)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl
                << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }
  return 0;

}