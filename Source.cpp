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

			rev -.05v
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
//dev generated
#include "fNum.h"//short hangs while generating large amounts of numbers, high CPU usage per thread.
using namespace std;
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
//####CLASS NEURON #######!!!!!!!!!!!!!!#######################|
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
struct Connection
{
	double weight;
	double deltaWeight;
};
//
class Neuron 
{ //if it's const it dosnt modify the obj
	public:
		Neuron(unsigned numOutputs, unsigned myIndex);
		void   setOutputVal (double val) { m_outputVal = val;  }
		double getOutputVal (void )const { return m_outputVal; }
		void   feedForward  ( const Layer& prevLayer );
		void calcOutputGradients(double targetVal);
		void calcHiddenGradients(const Layer &nextLayer); 
		void updateInputWeights(Layer& prevLayer);
		m_recentAverageError
		m_recentAverageErrorSmoothingFactor
			
	private:
		static double transferFunction(double x);
		static double transferFunctionDerivative(double x); //HYPERBOLIX TAN FUNC
		static double randomWeight(void) { return fRand(); }
		double sunDOW(const Layer& nextLayer);
		double m_outputVal;
		vector<Connection> m_outputWeights;
		unsigned m_myIndex;
		double   m_gradient;

};

void updateInputWeights(Layer& prevLayer) 
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
		Neuron& neuron = prevLayer[n];
		double oldDeltaWeight = neruon.m_outputWeights[m_myIndex].deltaWeight;
		double newDeltaWeight =
			// Individiual input, magnified by the gradient and train rate:
			eta
			* neuron_getOutputVal()
			* m_gradient
			// add momemntum = a fraction of the prev delta weight (alpha)
			+ alpha
			* oldDeltaWeight;


	}

}

double Neuron::sunDOW(const Layer& nextLayer)
{
	double sunDOW = 0.0;
	// sum contributions of the erros at the node we feed
	for (unsigned n=0; n < nextLayer.size() - 1; n++)
	{
		sum += m_outputWeights[n].weight * nextLayer[n].m_gradient;
	}
	return sum;
}

void Neuron::calcHiddenGradients(const Layer& nextLayer)
{ //figure out an error delta by the sum of the diravatives
  //of the weights of the next layer
	double dow = sunDOW(nextLayer);
	m_gradient = dow * Neuron::transferFunctionDerivative(m_outputVal);
	
}

void Neuron::calcOutputGradients(double targetVal);
{	//This keeps the training to where it reduces the overall net area
	double delta = targetVal - m_outputVal;
	m_gradient = delta * Neuron::transferFunctionsderiviative(m_outputVal);
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

Neuron::Neuron(unsigned numOutputs) //neruon consatructor
{	//c is for connections
	for (unsigned c = 0; c < numOutputs; c++)
	{
		//append (pushback) each itieration, new connection to the output weights container
		m_outputWeights.push_back(Connection());
		m_outputWeights.back().weight = randomWeight();
	} 
}

void Neuron::feedForward(const Layer& prevLayer) 
{
	double sum = 0.0;
	// Sum the prev layers output (whicher are the inputs)
	//Include bias node from prev layer

	for (unsigned n=0; n <prevLayer.size(); n++)
	{	//all inputs time(*) their connection weights
		sum += prevLayer[n]getOutputVal() *
			   prevLayer[n].m_outputWeights[m_myIndex].weight;
	}
	m_outputVal = Neuron::transferFunction(sum);
}
typedef vector<Neuron> Layer;

//------------------------------------------------------------|
//****NET CLASS DECLARTION************************************|
//------------------------------------------------------------|
class Net 
{
	public:
		Net(const vector<unsigned> &topology);
		void feedForward (const vector<double> &inputVals ); //get 'inputVal by refrence to a vector of doubles instead of copy and without changing anything
		void backProp    (const vector<double> &targetVals); //this is just a happy little sammich isnt it?
		void getResults  (		vector<double> &resultVals) const; //note the location of const here! 
	private:
		vector<Layer> m_layers; // m_layers[layerNum][neruonNum] a private data member
		double        m_error ; //this is RMS!!!
};
//  |~ CREATING NET
//  |
/***V*********************************************************/
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
		for (unsigned neruonNum = 0; neruonNum <= topology[layerNum]; ++neruonNum)
		{
		 //add one extra bias neuron into additon to neurons (<=)
			m_layers.back().push_back(Neuron(numOutputs));
			cout << "Neuron Created." << endl;
		}
	}	
}
//...........
void Net::feedForward(const vector<double>& inputVals)
{
	assert(inputVals.size() == m_layers.size() - 1);
	//assign (latach) the input values into the input neruons
	for (unsigned i = 0; i < inputVals.size(); ++i) {
		m_layers[0][i].setOutputVal(inputVals[i])); // 'i'th element
	}

	//forward prep
	for (unsigned layerNum = 1; layerNum < m_layers.size(); ++layerNum) {
		Layer& prevLayer = m_layers[layerNum - 1];
		for (unsigned n = 0; n < m_layers[layerNum].size() - 1; n++) {
			m_layers[layerNum][n].feedForward(prevLayer);
		}
	}
}
//
void Net::backProp(const vector<double> &targetVals)
{/* Calculate 
			|- net error. (RMS)
            |- output layer gradients.
			|- gradients on hidden layers.
			|- For every layer, from output to 1st hidden layer
			|- & update connection weights.*/
	Layer &outputLayer = m_layers.back();
	m_error = 0.0;
	for (unsigned n = 0; n < outputLayer.size() - 1; n++)
	{
		double delta = targetVals[n] - outputLayer[n].getOutputVal();
		m_error += delta * delta;
	}
	m_error /= outputLayer.size() - 1; //get average error ^2
	m_error = sqrt(m_error);           //RootMeanSquared
	//...................................
	//Implement a recent avg measurement:
	m_recentAverageError =
		(m_recentAverageError * m_recentAverageErrorSmoothingFacotr + m_error) /
		(m_recentAverageErrorSmoothingFactor + 1.0);
	//output gradients
	for (unsigned n = 0; N < outputLayer.size() - 1; n++)
	{
		outputLayer[n].calcOuputGradients(targetVals[n]);
	}
	//calc output gradients (on hidden layers)
	for (unsigned layerNum  = m_layers.size() - 2; layerNum > 0; layerNum--)
	{
		Layer& hiddenLayer = m_layers[layerNum    ];
		Layer& nextLayer   = m_layers[layerNum + 1];
		//for each one well add...
		for (unsigned n = 0; n < hiddenLayer.size(); n++)
		{
			hiddenLayer[n].calcHiddenGradients(nextLayer);
		}

	}
	//For all layers from output to first hidden layer,
	//		update connections weights.

	for (unsigned layerNum = m_layers.size() - 1; layerNum > 0; layerNum)
	{
		layer[n].updateInputWeights(prevLayer);
	}
}
//
void Net::getResults(vector<double>& resultVals) const
{
}
/*******end of net.*****************************************/

int main()
{

	while (true)
	{
		
		cout << fRand() << endl;
		
	}

	//{ 3   ,   2   ,   1} example
	vector <unsigned> topology;
	topology.push_back(3); //this is bad practicce (repeating code)
	topology.push_back(2); // 3 - 2 -1 = 3 inputs, 1 output
	topology.push_back(1); // 1 bias neuron per layer ( 3 layers)
	Net myNet(topology);


	

/*************************************************************/
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
	return 0;
}



