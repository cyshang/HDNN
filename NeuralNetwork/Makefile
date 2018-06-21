objects = main.o Parameter.o NeuralNetwork.o NetworkInfo.o Optimizer.o PES.o SDBP.o LMBP.o EKF.o \
          GroupBase.o GroupSDBP.o GroupEKF.o GroupLMBP.o LayerBase.o InputLayer.o HiddenLayer.o OutputLayer.o \
		  SymFunction.o FunctionInfo.o Molecule.o MonteCarloSetting.o Atom.o FuncType.o
EIGEN = -I ../Tools/Eigen3/
LIB = -L /work1/soft/intel2015/composer_xe_2015.2.164/mkl/lib/intel64 -lmkl_rt
MKL = -I /work1/soft/intel2015/mkl/include
options = -O3

TARGET = NeuralNetwork

$(TARGET) : $(objects)
	g++ -o SymFunction $(objects) -lm $(LIB)

main.o : main.cpp global.h NeuralNetwork.h SymFunction.h
	g++ $(EIGEN) $(MKL) -c main.cpp $(options)

Parameter.o : Parameter.cpp Parameter.h
	g++ $(EIGEN) $(MKL) -c Parameter.cpp $(options)	
	
NeuralNetwork.o : NeuralNetwork.cpp NeuralNetwork.h NetworkInfo.h GroupBase.h LayerBase.h Optimizer.h PES.h SDBP.h LMBP.h EKF.h
	g++ $(EIGEN) $(MKL) -c NeuralNetwork.cpp $(options)

NetworkInfo.o : NetworkInfo.cpp NetworkInfo.h
	g++ $(EIGEN) $(MKL) -c NetworkInfo.cpp $(options)

Optimizer.o : Optimizer.cpp Optimizer.h NeuralNetwork.h NetworkInfo.h GroupBase.h
	g++ $(EIGEN) $(MKL) -c Optimizer.cpp $(options)

PES.o : PES.cpp PES.h GroupBase.h NeuralNetwork.h NetworkInfo.h
	g++ $(EIGEN) $(MKL) -c PES.cpp $(options)

SDBP.o : SDBP.cpp SDBP.h NeuralNetwork.h NetworkInfo.h GroupSDBP.h
	g++ $(EIGEN) $(MKL) -c SDBP.cpp $(options)

LMBP.o : LMBP.cpp LMBP.h NeuralNetwork.h NetworkInfo.h GroupLMBP.h
	g++ $(EIGEN) $(MKL) -c LMBP.cpp $(options)

EKF.o : EKF.cpp EKF.h NeuralNetwork.h NetworkInfo.h GroupEKF.h
	g++ $(EIGEN) $(MKL) -c EKF.cpp $(options)

GroupBase.o : GroupBase.cpp GroupBase.h NeuralNetwork.h NetworkInfo.h
	g++ $(EIGEN) $(MKL) -c GroupBase.cpp $(options)

GroupSDBP.o : GroupSDBP.cpp GroupSDBP.h SDBP.h NeuralNetwork.h
	g++ $(EIGEN) $(MKL) -c GroupSDBP.cpp $(options)	
	
GroupEKF.o : GroupEKF.cpp GroupEKF.h EKF.h NeuralNetwork.h
	g++ $(EIGEN) $(MKL) -c GroupEKF.cpp $(options)

GroupLMBP.o : GroupLMBP.cpp GroupLMBP.h LMBP.h NeuralNetwork.h
	g++ $(EIGEN) $(MKL) -c GroupLMBP.cpp $(options)

LayerBase.o : LayerBase.cpp LayerBase.h NetworkInfo.h GroupBase.h
	g++ $(EIGEN) $(MKL) -c LayerBase.cpp $(options)
	
InputLayer.o : InputLayer.cpp InputLayer.h GroupBase.h
	g++ $(EIGEN) $(MKL) -c InputLayer.cpp $(options)

HiddenLayer.o : HiddenLayer.cpp HiddenLayer.h GroupBase.h
	g++ $(EIGEN) $(MKL) -c HiddenLayer.cpp $(options)

OutputLayer.o : OutputLayer.cpp OutputLayer.h GroupBase.h
	g++ $(EIGEN) $(MKL) -c OutputLayer.cpp $(options)

SymFunction.o : SymFunction.cpp SymFunction.h FunctionInfo.h Molecule.h FuncType.h NeuralNetwork.h MonteCarloSetting.h
	g++ $(EIGEN) $(MKL) -c SymFunction.cpp $(options)

FunctionInfo.o : FunctionInfo.cpp FunctionInfo.h
	g++ $(EIGEN) $(MKL) -c FunctionInfo.cpp $(options)

Molecule.o : Molecule.cpp Molecule.h
	g++ $(EIGEN) $(MKL) -c Molecule.cpp $(options)
	
MonteCarloSetting.o : MonteCarloSetting.cpp MonteCarloSetting.h
	g++ $(EIGEN) $(MKL) -c MonteCarloSetting.cpp $(options)

Atom.o : Atom.cpp Atom.h Molecule.h FunctionInfo.h SymFunction.h
	g++ $(EIGEN) $(MKL) -c Atom.cpp $(options)

FuncType.o : FuncType.cpp FuncType.h MonteCarloSetting.h
	g++ $(EIGEN) $(MKL) -c FuncType.cpp $(options)
	
.PHONY : clean
clean :
	-rm $(TARGET) $(objects)
