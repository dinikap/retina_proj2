#ifndef RETINA_PROJ2_H_
#define RETINA_PROJ2_H_

#include <vector>
#include <math.h>
#include "biodynamo.h"
#include "../../biodynamo/src/core/substance_initializers.h"
//#include "math_util.h"
//#include "substance_initializers.h"

/*
MODIFIED - to allow this to run in the latest BDM
Dinika P.
*/

namespace bdm {

/* Extend into cell class to define the cell types

MODIFIED
*/
  BDM_SIM_OBJECT(MyCell, Cell) {
    //BDM_SIM_OBJECT_HEADER(MyCellExt, 1, cellType);
    BDM_SIM_OBJECT_HEADER(MyCell, Cell, 1, cellType);

  public:
    MyCellExt() {}
    explicit MyCellExt(const std::array<double, 3>& position) : Base(position) {}

    //getter and setter for new data member
    void SetCellType(int t) { cellType[kIdx] = t; }
    int GetCellType() { return cellType[kIdx]; }
    // This function is used by ParaView for coloring the cells by their type
    int* GetCellTypePtr() { return cellType.data(); }

  private:
    //declare  new data member
    vec<int> cellType;
  };

  /*
  External substance will lead cell migration based on the concentration
  levels. Substance behavior is stated in a biology module
  */
 enum Substances { kSubstance };

  /*
    Created own function to create cells
    so all cells will start from a fixed bound and then migrate
    based on the substance concentration

    MODIFIED to work with latest version of BDM
    ref: model_initialiser.h, CreateCellsRandom
  */

 //template <typename Function, typename TResourceManager = ResourceManager<>>
 template <typename Function, typename TSimulation = Simulation<>>
 static void MyCellCreator(double min, double max, int num_cells, Function cell_builder) {
   auto* sim = TSimulation::GetActive();
   auto* rm = sim->GetResourceManager();
   auto* random = sim->GetRandom();
   // Determine simulation object type which is returned by the cell_builder
   using FunctionReturnType = decltype(cell_builder({0, 0, 0}));

   rm->template Reserve<FunctionReturnType>(num_cells);
   // so cells will be created at random only on the x and y axis
   // z axis is used to move cells to final resting position
   for (int i = 0; i < num_cells; i++) {
     double x = random->Uniform(min, max);
     double y = random->Uniform(min, max);
     //stop cells from moving in the z axis when generated
     double z = 0;
     auto new_simulation_object = cell_builder({x, y, z});
     rm->push_back(new_simulation_object);
   }
 }


/*
  BiologyModules for each cell type - specifies the behaviours of each cell type
  MODIFIED - the way BiologyModules are declared has changed
*/

/*
  Ganglion layer cells -> Assume only 1 type for the simplicity of model
  Ganglion cells will link to amacrine and bipolar cells
*/
    struct ganglionCell : public BaseBiologyModule {
      ganglionCell() : BaseBiologyModule(gAllEventIds){}

    //  template <typename T>
    template <typename TEvent, typename TBm>
    ganglionCell(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

    //add this line to declare T
    template <typename T, typename TSimulation = Simulation<>>
      void Run(T* cell){
        //initialisation of diffusion grid
        // if (!init_) {
        //   dg_= GetDiffusionGrid(kSubstance);
        //   init_ = true;
        // }
        //Edit
        auto* sim = TSimulation::GetActive();
        auto* rm = sim->GetResourceManager();
        dg_ = rm->GetDiffusionGrid(kSubstance);

        //get the position of cell and grad
        auto& position = cell->GetPosition();
        dg_->GetGradient(position, &gradient_);
        double concentration = dg_->GetConcentration(position);
        //use movement to make the cells move slower
        movement[0] = gradient_[0]*0.5;
        movement[1] = gradient_[1]*0.5;
        movement[2] = gradient_[2]*0.5;

        if(concentration < 0.0000025){
          cell->UpdatePosition(movement);
          //cell->SetPosition(cell->GetMassLocation());
          cell->SetPosition(cell->GetPosition());
        }
      }
      private:
        //bool init_ = false;
        DiffusionGrid* dg_ = nullptr;
        std::array<double, 3> gradient_;
        std::array<double, 3> movement;

        //BDM_CLASS_DEF_NV(ganglionCell, 1);
    };

    /*
    Amacrine cells that are in the INNER LIMITING LAYER
    link bipolar and ganglion cells
    cell takes input from ganglion cell to bipolar cell
    */
    struct amacrineCell : public BaseBiologyModule {
      amacrineCell() : BaseBiologyModule(gAllEventIds){}

      //template <typename T>
      template <typename TEvent, typename TBm>
      amacrineCell(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

      //add this line to declare T
      template <typename T, typename TSimulation = Simulation<>>
      void Run(T* cell){
        //initialisation of diffusion grid
        // if (!init_) {
        //   dg_= GetDiffusionGrid(kSubstance);
        //   init_ = true;
        // }
        //Edit
        auto* sim = TSimulation::GetActive();
        auto* rm = sim->GetResourceManager();
        dg_ = rm->GetDiffusionGrid(kSubstance);

        //get the position of cell and grad
        auto& position = cell->GetPosition();
        dg_->GetGradient(position, &gradient_);
        double concentration = dg_->GetConcentration(position);
        //use movement to make the cells move slower
        movement[0] = gradient_[0]*0.5;
        movement[1] = gradient_[1]*0.5;
        movement[2] = gradient_[2]*0.5;

          if(concentration < 0.0000002){
          cell->UpdatePosition(movement);
          //cell->SetPosition(cell->GetMassLocation());
          cell->SetPosition(cell->GetPosition());
        }
      }
    private:
      //bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
      std::array<double, 3> movement;
    };

    /*
    Bipolar cells
    exists between photoreceptors and ganglion cells
    so it will either synpase with:
    1. photoreceptors -> rods/ cones (i.e. Parasol cells)
    2. horizontal cells
    */
    struct bipolarCell : public BaseBiologyModule {
      bipolarCell() : BaseBiologyModule(gAllEventIds){}

      //template <typename T>
      template <typename TEvent, typename TBm>
      bipolarCell(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

      //add this line to declare T
      template <typename T, typename TSimulation = Simulation<>>
      void Run(T* cell){
        //initialisation of diffusion grid
        // if (!init_) {
        //   dg_= GetDiffusionGrid(kSubstance);
        //   init_ = true;
        // }
        //Edit
        auto* sim = TSimulation::GetActive();
        auto* rm = sim->GetResourceManager();
        dg_ = rm->GetDiffusionGrid(kSubstance);

        //get the position of cell and grad
        auto& position = cell->GetPosition();
        dg_->GetGradient(position, &gradient_);
        double concentration = dg_->GetConcentration(position);
        //use movement to make the cells move slower
        movement[0] = gradient_[0]*0.5;
        movement[1] = gradient_[1]*0.5;
        movement[2] = gradient_[2]*0.5;

        if (concentration < 0.000000025) {
          cell->UpdatePosition(movement);
          //cell->SetPosition(cell->GetMassLocation());
          cell->SetPosition(cell->GetPosition());
        }
      }
    private:
      //bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
      std::array<double, 3> movement;
    };

    /*horizontal cells
    part 2:
    cells work laterally
    connected to outputs from rods and cones
    */
    struct horizontalCell : public BaseBiologyModule {
      horizontalCell() : BaseBiologyModule(gAllEventIds){}

      //template <typename T>
      template <typename TEvent, typename TBm>
      horizontalCell(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

      //add this line to declare T
      template <typename T, typename TSimulation = Simulation<>>
      void Run(T* cell){
        //initialisation of diffusion grid
        // if (!init_) {
        //   dg_= GetDiffusionGrid(kSubstance);
        //   init_ = true;
        // }
        //Edit
        auto* sim = TSimulation::GetActive();
        auto* rm = sim->GetResourceManager();
        dg_ = rm->GetDiffusionGrid(kSubstance);

        //get the position of cell and grad
        auto& position = cell->GetPosition();
        dg_->GetGradient(position, &gradient_);
        double concentration = dg_->GetConcentration(position);
        //use movement to make the cells move slower
        movement[0] = gradient_[0]*0.5;
        movement[1] = gradient_[1]*0.5;
        movement[2] = gradient_[2]*0.5;

        //if (concentration < 0.000000037) {
        if (concentration < 0.000000007) {
          cell->UpdatePosition(movement);
          //cell->SetPosition(cell->GetMassLocation());
          cell->SetPosition(cell->GetPosition());
        }
      }
    private:
      //bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
      std::array<double, 3> movement;
    };

    /*cones
    part 2:
    cells are long but rounder and wider than rods
    connect to horizontal cells + bipolar cells
    */
    struct coneCell : public BaseBiologyModule {
      coneCell() : BaseBiologyModule(gAllEventIds){}

      //template <typename T>
      template <typename TEvent, typename TBm>
      coneCell(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

      //add this line to declare T
      template <typename T, typename TSimulation = Simulation<>>
      void Run(T* cell){
        //initialisation of diffusion grid
        // if (!init_) {
        //   dg_= GetDiffusionGrid(kSubstance);
        //   init_ = true;
        // }
        //Edit
        auto* sim = TSimulation::GetActive();
        auto* rm = sim->GetResourceManager();
        dg_ = rm->GetDiffusionGrid(kSubstance);

        //get the position of cell and grad
        auto& position = cell->GetPosition();
        dg_->GetGradient(position, &gradient_);
        double concentration = dg_->GetConcentration(position);
        //use movement to make the cells move slower
        movement[0] = gradient_[0]*0.5;
        movement[1] = gradient_[1]*0.5;
        movement[2] = gradient_[2]*0.5;

        //cones do not migrate far off from basal
        if (concentration < 0.000000003) {
          cell->UpdatePosition(movement);
        //  cell->SetPosition(cell->GetMassLocation());
        cell->SetPosition(cell->GetPosition());
        }
      }
    private:
      //bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
      std::array<double, 3> movement;
    };

    /*rods
    part 2:
    cells are longish -> outer seg + inner seg + nucleus
    connect to horizontal cells + bipolar cells
    */
    struct rodCell : public BaseBiologyModule {
      rodCell() : BaseBiologyModule(gAllEventIds){}

      //template <typename T>
      template <typename TEvent, typename TBm>
      rodCell(const TEvent& event, TBm* other, uint64_t new_oid = 0) {}

      //add this line to declare T
      template <typename T, typename TSimulation = Simulation<>>
      void Run(T* cell){
        //initialisation of diffusion grid
        // if (!init_) {
        //   dg_= GetDiffusionGrid(kSubstance);
        //   init_ = true;
        // }
        //Edit
        auto* sim = TSimulation::GetActive();
        auto* rm = sim->GetResourceManager();
        dg_ = rm->GetDiffusionGrid(kSubstance);

        //get the position of cell and grad
        auto& position = cell->GetPosition();
        dg_->GetGradient(position, &gradient_);
        double concentration = dg_->GetConcentration(position);
        //use movement to make the cells move slower
        movement[0] = gradient_[0]*0.5;
        movement[1] = gradient_[1]*0.5;
      //  movement[2] = gradient_[2]*0.5;

        //rods do not migrate far off from basal
        if (concentration < 0.000000003) {
          cell->UpdatePosition(movement);
        //  cell->SetPosition(cell->GetMassLocation());
        cell->SetPosition(cell->GetPosition());
        }
      }
    private:
      //bool init_ = false;
      DiffusionGrid* dg_ = nullptr;
      std::array<double, 3> gradient_;
      std::array<double, 3> movement;
    };

// Define compile time parameter
// template <typename Backend>
// struct CompileTimeParam : public DefaultCompileTimeParam<Backend> {
//   using BiologyModules = Variant<ganglionCell, amacrineCell, bipolarCell, horizontalCell, coneCell, rodCell>;
//   using AtomicTypes = VariadicTypedef <MyCell>;
// };

  BDM_CTPARAM() {
    BDM_CTPARAM_HEADER();
    using SimObjectTypes = CTList<MyCell>;

    //override default biology modules for cell
    BDM_CTPARAM_FOR(bdm, MyCell) { using BiologyModules = CTList<ganglionCell, amacrineCell, bipolarCell, horizontalCell, coneCell, rodCell>;};
  };

inline int Simulate(int argc, const char** argv) {
  //InitializeBioDynamo(argc, argv);

  /*
    Params for Paraview
  */
  // Param::bound_space_ = true;
  // Param::min_bound_ = 0;
  // //max bound is 250um
  // Param::max_bound_ = 250;
  // Param::run_mechanical_interactions_ = true;
  auto set_param = [] (auto* param) {
    param->bound_space_ = true;
    param->min_bound_ = 0;
    param->max_bound_ = 250;
    param->run_mechanical_interactions_ = true;
    param->live_visualization_ = true;
    param->export_visualization_ = true;
    param->visualization_export_interval_ = 2;
    param->visualize_sim_objects_["MyCell"] = std::set<std::string>{"cellType"};
  };

   Simulation<> simulation(argc, argv, set_param);
  // auto* rm = simulation.GetResourceManager();
   auto* param = simulation.GetParam();
   auto* random = simulation.GetRandom();
  /* this is to ensure that the model is reproducible for a
  specific seed
  */
//  int randSeed = rand() % 1000;    gTRandom.SetSeed(randSeed);
  cout << "modelling seed: " << random <<endl;

  /* Define cells to simulate:
  the cell diameters, cell type and name will be defined here
  */

  auto construct_ganglion = [](const std::array<double, 3>& position) {
        MyCell cell(position);
        //estimate gangliong cell diameter to be 11um
        cell.SetDiameter(11);
        cell.AddBiologyModule(ganglionCell());
        cell.SetCellType(6);
        return cell;
      };
    cout << "Ganglion cells created" << endl;
    MyCellCreator(param->min_bound_, param->max_bound_, 400, construct_ganglion);

    auto construct_amacrine = [](const std::array<double, 3>& position){
        MyCell cell(position);
        //assume average dendritic field to be small/medium field as we are looking
        //at central ret so around 100um
        //170418 --> ignore the dendtritic fields, estimate cell diameter
        //approx 8.75
        cell.SetDiameter(9);
        cell.AddBiologyModule(amacrineCell());
        cell.SetCellType(5);
        return cell;
      };
    cout << "Amacrine cells created" << endl;
    MyCellCreator(param->min_bound_, param->max_bound_, 400, construct_amacrine);

      auto construct_bipolar = [](const std::array<double, 3>& position){
        MyCell cell(position);
        //assume avg diamter using midget cone bipolar fmB and imB of 12um in 4.5mm region

        cell.SetDiameter(9);
        cell.AddBiologyModule(bipolarCell());
        cell.SetCellType(4);
        return cell;
      };
      cout << "Bipolar cells created" << endl;
      MyCellCreator(param->min_bound_, param->max_bound_, 400, construct_bipolar);

      auto construct_horizontal = [](const std::array<double, 3>& position){
          MyCell cell(position);
          //assume avg diamter using H1 of 25um dentritic field in 2.5mm region
          //170418 --> changed to estimated cell diameter 7.5
          cell.SetDiameter(8);
          cell.AddBiologyModule(horizontalCell());
          cell.SetCellType(3);
          return cell;
      };
      cout << "Horizontal cells created" << endl;
      MyCellCreator(param->min_bound_, param->max_bound_, 200, construct_horizontal);

      auto construct_cone = [](const std::array<double, 3>& position){
        MyCell cell(position);
        //approximately 2um in foveal area
        cell.SetDiameter(2);
        cell.AddBiologyModule(coneCell());
        cell.SetCellType(1);
        return cell;
      };
      cout << "Cone cells created" << endl;
      MyCellCreator(param->min_bound_, param->max_bound_, 250, construct_cone);

      auto construct_rod = [](const std::array<double, 3>& position){
        MyCell cell(position);
        //approximately 2umn in diameter
        cell.SetDiameter(2);
        cell.AddBiologyModule(rodCell());
        cell.SetCellType(2);
        return cell;
      };
      cout << "Rod cells created" << endl;
      MyCellCreator(param->min_bound_, param->max_bound_, 250, construct_rod);
    //defining substances in simulation
    //diffusion coefficient of 0.5, a decay constant 0f 0.1 and a resolution of 1
    ModelInitializer::DefineSubstance(kSubstance, "kSubstance", 0.5, 0.1, 4);
    //initialise substance: enum of substance, name, function type used
    //mean value of 200 along the z-axis, and a variance of 100
    ModelInitializer::InitializeSubstance(kSubstance, GaussianBand(200, 100, Axis::kZAxis));


  //link to paraview to show visualization
  //EDIT: set this in set_param
    // Param::live_visualization_ = true;
    // Param::export_visualization_ = true;
    // Param::visualization_export_interval_ = 2;
    // Param::visualize_sim_objects_["MyCell"] = std::set<std::string>{"cellType"};

  // Run simulation for one timestep
  Scheduler<> scheduler;
  int maxStep = 1800;
  for (int i = 0; i < maxStep; i++){
    scheduler.Simulate(1);
  }


  std::cout << "Simulation completed successfully!" << std::endl;
  return 0;
}

} // namespace bdm

#endif // RETINA_PROJ2_H_
