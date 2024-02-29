// FUNCTION DECLARATION

// TGA FUNCTIONs
void OdeTGA(const double* y, const double t, double* dy, void* args);
void MyStepPrintTGA(const double t, const Eigen::VectorXd& x);
void PrintFinalStatusTGA(const Eigen::VectorXd& final_status);

// 1D FUNCTIONs
void OdeTotal(const double* y, const double t, double* dy, void* args);
void MyStepPrintTotal(const double t, const Eigen::VectorXd& x);
void PrintFinalStatusTotal(const Eigen::VectorXd& final_status);

void ChangeDimensionsFunction();
void CreateFractions();
void SetVolume(const double t);
void CalculateProperties();
void CalculateKinetics();
void CreateMaterialBalances();
void CreateEnergyBalances();
