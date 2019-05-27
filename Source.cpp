#pragma  once
/*
origin is azure
master is gitHub
	ITJ
	may 21 2019
	NeuralNet Prototype
	https://www.youtube.com/watch?v=KkwX7FkLfug

	"TEST EARLY.
		   TEST, OFTEN."
			

			GOAL: Designa  function Deep Neural net and understand basic concepts. Thenb make a way too feed it constant data
			after this project is complete I will then port this to my FPGA in RTL

			Trident layers pays after all!

			rev version  1.1a
*/
//#############################################################|


#include <iostream>
#include <cmath>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <chrono>
#include <random>
#include <cassert>
#include <fstream>
#include <sstream>
#include <cstdlib>
//dev generated
#include "fNum.h"//short hangs while generating large amounts of numbers, high CPU usage per thread.
using namespace std;
//#############################################################|
//~~~x~~x~~xtemp class to figure stuff out with~x~~x~~~x~~~x~~cp
//#############################################################|

class TrainingData
{
public:
	TrainingData(const string filename);
	bool isEof(void) { return m_trainingDataFile.eof(); }
	void getTopology(vector<unsigned>& topology);

	// Returns the number of input values read from the file:
	unsigned getNextInputs(vector<double>& inputVals);
	unsigned getTargetOutputs(vector<double>& targetOutputVals);

private:
	ifstream m_trainingDataFile;
};

void TrainingData::getTopology(vector<unsigned>& topology)
{
	string line;
	string label;

	getline(m_trainingDataFile, line);
	stringstream ss(line);
	ss >> label;
	if (this->isEof() || label.compare("topology:") != 0) {
		abort();
	}

	while (!ss.eof()) {
		unsigned n;
		ss >> n;
		topology.push_back(n);
	}

	return;
}

TrainingData::TrainingData(const string filename)
{
	m_trainingDataFile.open(filename.c_str());
}

unsigned TrainingData::getNextInputs(vector<double>& inputVals)
{
	inputVals.clear();

	string line;
	getline(m_trainingDataFile, line);
	stringstream ss(line);

	string label;
	ss >> label;
	if (label.compare("in:") == 0) {
		double oneValue;
		while (ss >> oneValue) {
			inputVals.push_back(oneValue);
		}
	}

	return inputVals.size();
}

unsigned TrainingData::getTargetOutputs(vector<double>& targetOutputVals)
{
	targetOutputVals.clear();

	string line;
	getline(m_trainingDataFile, line);
	stringstream ss(line);

	string label;
	ss >> label;
	if (label.compare("out:") == 0) {
		double oneValue;
		while (ss >> oneValue) {
			targetOutputVals.push_back(oneValue);
		}
	}

	return targetOutputVals.size();
}


//~~~~~~x~~x~~~~x~~~~~x~~~x~~~x~~~x~~~~x~~~~~x~~~x~~x~~~~xx~~~xcp

//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//####CLASS NEURON #######!!!!!!!!!!!!!!#######################|
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
struct Connection
{
	double weight;
	double deltaWeight;
};
//..............
class Neuron; //THIS LINE IS SUPER IMPORTANT (duh)
typedef vector<Neuron> Layer;


class Neuron 
{ //if it's const it dosnt modify the obj
public:
	Neuron(unsigned numOutputs, unsigned myIndex);
	void setOutputVal(double val) { m_outputVal = val; }
	double getOutputVal(void) const { return m_outputVal; }
	void feedForward(const Layer& prevLayer);
	void calcOutputGradients(double targetVal);
	void calcHiddenGradients(const Layer& nextLayer);
	void updateInputWeights(Layer& prevLayer);
private:
	static double eta;   // [0.0..1.0] overall net training rate
	static double alpha; // [0.0..n] multiplier of last weight change (momentum)
	static double transferFunction(double x);
	static double transferFunctionDerivative(double x);
		static double randomWeight(void) { return fRand(); }
		double sumDOW(const Layer& nextLayer) const;
		double m_outputVal;
		vector<Connection> m_outputWeights;
		unsigned m_myIndex;
		double m_gradient;
};

double Neuron::eta = 0.15;    // overall net learning rate, [0.0..1.0]
double Neuron::alpha = 0.5;   // momentum, multiplier of last deltaWeight, [0.0..1.0]


void Neuron::updateInputWeights(Layer& prevLayer)
{
	// the weights to be updated are in the connection container
	// in the neruons in the preceding layer
	/*
		(eta)   - overall net learning rate
				   0.0 - slow     learner
				   0.2 - medium   learner
				   1.0 - reckless learner
		(alpha) - momentum
				   0.0 - NO       momentum
				   0.5 - moderate momentum
				   1.0 - broken   momentum
	*/


	for (unsigned n = 0; n < prevLayer.size(); n++)
	{
		for (unsigned n = 0; n < prevLayer.size(); ++n) {
			Neuron& neuron = prevLayer[n];
			double oldDeltaWeight = neuron.m_outputWeights[m_myIndex].deltaWeight;
			//.....................
			double newDeltaWeight =
				// Individiual input, magnified by the gradient and train rate:
				eta
				* neuron.getOutputVal()
				* m_gradient
				// add momemntum = a fraction of the prev delta weight (alpha)
				+ alpha
				* oldDeltaWeight;

			neuron.m_outputWeights[m_myIndex].deltaWeight = newDeltaWeight;
			neuron.m_outputWeights[m_myIndex].weight += newDeltaWeight;
		}
	}
}
double Neuron::sumDOW(const Layer & nextLayer) const
{
		double sum = 0.0;

		// Sum our contributions of the errors at the nodes we feed.

		for (unsigned n = 0; n < nextLayer.size() - 1; ++n) {
			sum += m_outputWeights[n].weight * nextLayer[n].m_gradient;
		}

		return sum;
}
void Neuron::calcHiddenGradients(const Layer& nextLayer)
{ //figure out an error delta by the sum of the diravatives
  //of the weights of the next layer
	double dow = sumDOW(nextLayer);
	m_gradient = dow * Neuron::transferFunctionDerivative(m_outputVal);
}

void Neuron::calcOutputGradients(double targetVal)
{	//This keeps the training to where it reduces the overall net area
	double delta = targetVal - m_outputVal;
	m_gradient = delta * Neuron::transferFunctionDerivative(m_outputVal);
}

double Neuron::transferFunction(double x)
{
	//tanh - output range [-1....1]
	return tanh(x);
}

double Neuron::transferFunctionDerivative(double x)
{
	// tanh derivative
	return 1.0 - x * x;
}


Neuron::Neuron(unsigned numOutputs, unsigned myIndex)//neruon consatructor
{	//c is for connections
	for (unsigned c = 0; c < numOutputs; c++)
	{
		//append (pushback) each itieration, new connection to the output weights container
		m_outputWeights.push_back(Connection(       ));
		m_outputWeights.back().weight = randomWeight();
	} 
	m_myIndex = myIndex;
}

void Neuron::feedForward(const Layer& prevLayer) 
{
	double sum = 0.0;
	// Sum the prev layers output (whicher are the inputs)
	//Include bias node from prev layer

	for (unsigned n=0; n <prevLayer.size(); n++)
	{	//all inputs time(*) their connection weights
		sum += prevLayer[n].getOutputVal() *// ( times )
			   prevLayer[n].m_outputWeights[m_myIndex].weight;
	}
	m_outputVal = Neuron::transferFunction(sum);
}


//------------------------------------------------------------|
//****NET CLASS DECLARTION************************************|
//------------------------------------------------------------|
class Net 
{
	public:
		Net(const vector<unsigned> &topology);
		void feedForward (const vector<double> &inputVals ); //get 'inputVal by refrence to a vector of doubles instead of copy and without changing anything
		void backProp(const vector<double>& targetVals);  //this is just a happy little sammich isnt it?
		void getResults  (		vector<double> &resultVals) const; //note the location of const here! 
		double getRecentAverageError(void) const { return m_recentAverageError; }
	private:
		vector<Layer> m_layers; // m_layers[layerNum][neruonNum] a private data member
		double        m_error ; //this is RMS!!!
		double m_recentAverageError;
		static double m_recentAverageSmoothingFactor;
};
//  |~ CREATING NET
//  |
/***V*********************************************************/
double Net::m_recentAverageSmoothingFactor = 100.0; // Number of training samples to average over4
void Net::getResults(vector<double>& resultVals) const
{
	resultVals.clear();

	for (unsigned n = 0; n < m_layers.back().size() - 1; ++n) {
		resultVals.push_back(m_layers.back()[n].getOutputVal());
	}
}



Net::Net(const vector<unsigned> &topology)
{
	unsigned numLayers = topology.size();
	for (unsigned layerNum = 0;  layerNum < numLayers; ++layerNum) 
	{
		m_layers.push_back(Layer()); 
		unsigned numOutputs = layerNum == topology.size() - 1 ? 0 : topology[layerNum + 1];
	 //push back std container member function to append an element on,
	 //	a standerd container. (layer)	
	 //Made a new layer, now fill it with (n)ith neruons, add A bias,
	 //neuron to the Layer.
		for (unsigned neruonNum = 0; neruonNum <= topology[layerNum]; neruonNum++)
		{
		 //add one extra bias neuron into additon to neurons (<=)
			m_layers.back().push_back(Neuron(numOutputs, neruonNum));
			std::cout << "Made a Neuron!" << endl;
		}
		// Force the bias node's output to 1.0 (it was the last neuron pushed in this layer):
		m_layers.back().back().setOutputVal(1.0);
	}	
}
//...........
void Net::feedForward(const vector<double>& inputVals)
{
	assert(inputVals.size() == m_layers.size() - 1);
	//assign (latach) the input values into the input neruons
	for (unsigned i = 0; i < inputVals.size(); ++i) {
		m_layers[0][i].setOutputVal(inputVals[i]); // 'i'th element
	}

	//forward prep
	// forward propagate
	for (unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum) {
		Layer& prevLayer = m_layers[layerNum - 1];
		for (unsigned n = 0; n < m_layers[layerNum].size() - 1; ++n) {
			m_layers[layerNum][n].feedForward(prevLayer);
		}
	}
}
//
void Net::backProp(const vector<double>& targetVals)
{/* Calculate
			|- net error. (RMS)
			|- output layer gradients.
			|- gradients on hidden layers.
			|- For every layer, from output to 1st hidden layer
			|- & update connection weights.*/
	Layer& outputLayer = m_layers.back();
	m_error = 0.0;
	//............
	for (unsigned n = 0; n < outputLayer.size() - 1; ++n) {
		double delta = targetVals[n] - outputLayer[n].getOutputVal();
		m_error += delta * delta;
	}
	m_error /= outputLayer.size() - 1;//get average error ^2
	m_error = sqrt(m_error);           //RootMeanSquared
	//...................................
	//Implement a recent avg measurement:
	m_recentAverageError =
		(m_recentAverageError * m_recentAverageSmoothingFactor + m_error)
		/ (m_recentAverageSmoothingFactor + 1.0);
	//output gradients
	for (unsigned n = 0; n < outputLayer.size() - 1; n++)
	{
		outputLayer[n].calcOutputGradients(targetVals[n]);
	}
	//calc output gradients (on hidden layers)
	for (unsigned layerNum = m_layers.size() - 2; layerNum > 0; layerNum--)
	{
		Layer& hiddenLayer = m_layers[layerNum];
		Layer& nextLayer = m_layers[layerNum + 1];
		//for each one well add...
		for (unsigned n = 0; n < hiddenLayer.size(); n++)
		{
			hiddenLayer[n].calcHiddenGradients(nextLayer);
		}
	}//moved a } from here to bottem of backprop
		
		//For all layers from output to first hidden layer,
		//		update connections weights.

		for (unsigned layerNum = m_layers.size() - 1; layerNum > 0; --layerNum) {
			Layer& layer = m_layers[layerNum];
			Layer& prevLayer = m_layers[layerNum - 1];

			for (unsigned n = 0; n < layer.size() - 1; ++n)
				{
					layer[n].updateInputWeights(prevLayer); //was 'layer'
				}//??? as was 'layer
	}
}

//

void showVectorVals(string label, vector<double>& v)
{
	cout << label << " ";
	for (unsigned i = 0; i < v.size(); ++i) {
		cout << v[i] << " ";
	}

	cout << endl;
}

/*******end of net.*****************************************/

int main()
{
	//test cp~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	TrainingData trainData("trainingData.txt");

	// e.g., { 3, 2, 1 }
	vector<unsigned> topology;
	trainData.getTopology(topology);

	Net myNet(topology);

	vector<double> inputVals, targetVals, resultVals;
	int trainingPass = 0;

	while (!trainData.isEof()) {
		++trainingPass;
		cout << endl << "Pass " << trainingPass;

		// Get new input data and feed it forward:
		if (trainData.getNextInputs(inputVals) != topology[0]) {
			break;
		}
		showVectorVals(": Inputs:", inputVals);
		myNet.feedForward(inputVals);

		// Collect the net's actual output results:
		myNet.getResults(resultVals);
		showVectorVals("Outputs:", resultVals);

		// Train the net what the outputs should have been:
		trainData.getTargetOutputs(targetVals);
		showVectorVals("Targets:", targetVals);
		assert(targetVals.size() == topology.back());

		myNet.backProp(targetVals);

		// Report how well the training is working, average over recent samples:
		cout << "Net recent average error: "
			<< myNet.getRecentAverageError() << endl;
	}

	cout << endl << "Done" << endl;
	//test cp~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	/*

	//{ 3   ,   2   ,   1} example
	vector <unsigned> topology;
	topology.push_back(3); //this is bad practicce (repeating code)
	topology.push_back(2); // 3 - 2 -1 = 3 inputs, 1 output
	topology.push_back(1); // 1 bias neuron per layer ( 3 layers)
	Net myNet(topology);


	

/*************************************************************
	//vector<double> topology;
	//--------------------------
	vector<double>    inputVals ;
	myNet.feedForward(inputVals); //appending?
	//...........................
	vector<double>   targetVals ;	
	myNet.backProp  (targetVals); //what the output should have been
	//...........................
	vector<double>   resultVals ;
	myNet.getResults(resultVals);
//|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
	*/
	return 0;
}



