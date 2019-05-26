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

			rev -.05v test
*/
//#############################################################|

#include <iostream>
#include <math.h>
#include <stdio.h>
#include <cstdio>
#include <vector>
#include <chrono>
#include <random>

#include "fNum.h"
//||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
using namespace std;
//####CLASS NEURON# ###########################################|
struct Connection
{
	double weight;
	double deltaWeight;
};
//
class Neuron 
{
	public:
		Neuron(unsigned numOutputs);
	private:
		static double randomWeight(void) { return fRand(); }
			double m_outputVal;
		vector<Connection> m_outputWeights;
		

};

Neuron::Neuron(unsigned numOutputs) //neruon consatructor
{	//c is for connections
	for (unsigned c = 0; c < numOutputs; c++)
	{
		//append (pushback) each itieration, new connection to the output weights container
		m_outputWeights.push_back(Connection());
		m_outputWeights.back().weight = randomWeight();
	} 
}

//------------------------------------------------------------|
typedef vector<Neuron> Layer;
/*****NET CLASS DECLARTION************************************/
class Net 
{
	public:
		Net(const vector<unsigned> &topology);
		void feedForward (const vector<double> &inputVals ); //get 'inputVal by refrence to a vector of doubles instead of copy and without changing anything
		void backProp    (const vector<double> &targetVals);
		void getResults  (		vector<double> &resultVals) const; //note the location of const here! 
	private:
		vector<Layer> m_layers; // m_layers[layerNum][neruonNum] a private data member
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
///////////////////////////////////////////////////////////////////
void Net::feedForward(const vector<double>& inputVals)
{
}
//
void Net::backProp(const vector<double>& backProp)
{
}
//
void Net::getResults(vector<double>& resultVals) const
{
}
/*************************************************************/
 float  sizer[2]; // testing var
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



