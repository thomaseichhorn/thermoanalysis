//C++ headers
#include <iostream>
#include <string>
#include <algorithm>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <stdio.h>

//Root headers
#include "TROOT.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TAxis.h"
#include "TMath.h"
#include "TLine.h"
#include "TLegend.h"
#include "TTree.h"
#include "TObject.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TFile.h"
#include "TDatime.h"
#include "TF1.h"
#include "TColor.h"
#include "TH1F.h"
#include "TPaveText.h"
#include "TMultiGraph.h"

// the namespace we are working in
using namespace std;

/*
Comments:
authors: Michael Bornholdt, Thomas Eichhorn

********************

Compile this script with:
g++ -I `root-config --incdir` -o test main.cc `root-config --libs` -Wall -std=c++0x -pedantic

Then run: ./test /path/to/runlist
Change test in the above two lines to a different executable name if you want.

********************

Run this script in ROOT with:

root -l -b
.x main.cc

Then input the runlist when prompted.

********************

Compile this script in ROOT with:

root -l -b
.x main.cc++

Then input the runlist when prompted.

********************

*/


// ********************
// runtime variables
// ********************

// the mode we want to run: 0 = testing, 1 = calibration, 2 = analysis, 3 = calibration and analysis
int mode = 3;

// output verbosity from loud = 0  to silent = 5
int debug = 3;

// the output file we want to save into
TFile * outputFile = new TFile("output.root", "RECREATE");

// write out every nth event
int precision = 50;

// the ntuple we read in
TTree *mytuple;

// the entry count of a tuple
long int tupleentrycount = 0;


// ********************
// constants:
// ********************

// number of sensors in the setup
const int sensors = 10;

// max allowed change in a sensor's temperature to be considered stable for calibration
const double deltacali = 0.001;

// max allowed change in a sensor's temperature to be considered stable for gradients
const double deltagrad = 0.0075;

// max number of sensor calibrations in a run
const int maxcalibs = 8;

// max number of stable temperature points for gradients
const int maxstable = 50;

// resistance of the heating element in ohms
const double resistor = 20.0;

// surface area of the contact between the aluminium blocks in m^2
// 30mm x 30mm
const double area = 0.000001 * 900;

// the effect of thermal grease, a temperature to be subtracted (twice) at the block interface
const double greasetemp = 0.00;

// assuming an error on all temperature measurements...
const double errorpercentage = 0.01;

// the total number of different measurement runs
const int maxmeas = 100;


// ********************
// the lists of files:
// ********************

// a list of the files to open
std::vector<std::string> filelist;

// the sorting of the sensors, from top to bottom
std::vector<std::string> sensorsort;

// list of broken sensors to skip
std::vector<std::string> brokenlist;

// the materials used in each measurement, sorted by measurement
std::vector<std::string> material;

// the different materials of all measurements
std::vector<std::string> materiallist;

// a count on the different materials
int materialcount;

// the thickness of the interface layer, sorted by measurement
std::vector<int> thickness;

// the different thicknesses of all measurements
std::vector<int> thicknesslist;

// a count on the different thicknesses
int thicknesscount;

// any other comments logged for a run
std::vector<std::string> comments;


// ********************
// the variables to read into:
// ********************

// the measurement time
unsigned int uTime;

// the sensor temperature
float temperature[sensors] = {0.0};

// the current applied to the resistor
float current1 = 0.0;

// the working temperature that is set
float workingTemperature = 0.0;


// ********************
// operational variables
// ********************

// ********************
// booleans
// ********************

// inside the calibration
bool insidecool, insideend, lookForCal = false;


// ********************
// floats
// ********************

// the working temperature at a calibration point
float work_Temperature[maxmeas][maxcalibs] = {0.0};

// the previous working temperature
float workingTemperaturebefore  = -100.0;

// the sensor temperature at the previous measurement point
float temperaturebefore[sensors] = {0.0};

// the resulting temperature difference in a sensor
float deltaT[sensors] = {0.0};

// the temperature difference between different sensors
float deltaDT[sensors] = {0.0};

// the temperature of a sensor at a calibration point
float calitemp[maxmeas][sensors][maxcalibs] = {0.0};

// the average temperature of a sensor over all calibrations
float avg_calitemp[maxmeas][sensors] = {0.0};

// points of stable temperature used for futher analysis
float stabletemp[maxmeas][sensors][maxstable] = {0.0};

// the working temperature at these points
float stablework[maxmeas][maxstable] = {0.0};

// the current at these points
float stablecurrent[maxmeas][maxstable] = {0.0};

// the time of these points
float stabletime[maxmeas][maxstable] = {0.0};

// the temperature difference between bottom and top blocks
float tempdiff[maxmeas][maxstable] = {0.0};

// the temperature at which the difference is measured
float tempdifftemp[maxmeas][maxstable] = {0.0};

// the average temperature difference
float avg_tempdiff = 0.0;

// the temperature gradient between both aluminium blocks
float gradient_blocks = 0.0;


// ********************
// integers
// ********************

// the number of points in a graph
int usedpoints = 0;

// the number of measurement points between two calibrations
int inbetween = 0;

// a temp time
int temptime = 0;

// the calibration counter
int calibs[maxmeas] = {0};

// the number of stable measurement points
int stablepoints[maxmeas] = {0};

// a sign to separate points in a plot
int pointsep = -1;


// ********************
// doubles
// ********************

// for the time
double time1, time0  = 0.0;

// the time of the calibration
double calitime[maxmeas][maxcalibs] = {0.0};


// ********************
// ROOTs
// ********************

// ********************
// calibration
// ********************

// graphs for the calibrations, also lines for the position in time
TGraphErrors* calibrationgraph[maxmeas][maxcalibs];
TLine* caliposition[maxmeas][maxcalibs];

// an average calibration graph
TGraphErrors* avg_calibrationgraph[maxmeas];

// a canvas for the calibrations
TCanvas* c_cali[maxmeas];

// a histogram for the calibrations
TH2D* h_cali[maxmeas];

// a legend for the calibrations
TLegend* l_cali[maxmeas];


// ********************
// temperature vs time
// ********************

// graphs for the sensor temperatures
TGraphErrors* tempgraph[maxmeas][sensors];
TGraphErrors* deltatempgraph[maxmeas][sensors];

// a canvas for the temperatures
TCanvas* c_temps[maxmeas];

// a histogram for the temperatures
TH2D* h_temps[maxmeas];

// a legend for the temperatures
TLegend* l_temps[maxmeas];


// ********************
// delta temperature vs time
// ********************

// a canvas for the delta temperatures
TCanvas* c_deltatemps[maxmeas];

// a histogram for the delta temperatures
TH2D* h_deltatemps[maxmeas];

// a legend for the delta temperatures
TLegend* l_deltatemps[maxmeas];


// ********************
// temperature gradients
// ********************

// graphs for the temperature gradients
// also fits for this
// and lines to show the position in time
TGraphErrors* gradgraph[maxmeas][maxstable];
TF1* gradfit1[maxmeas][maxstable];
TF1* gradfit2[maxmeas][maxstable];
TLine* gradposition[maxmeas][maxstable];

// a canvas for the temperature gradients
TCanvas* c_gradtemps[maxmeas];

// a histogram for the temperature gradients
TH2D* h_gradtemps[maxmeas];

// a legend for the temperature gradients
TLegend* l_gradtemps[maxmeas];

// a histogram for the temperature differences from the gradients
TH1D* h_blockdifference[maxmeas];


// ********************
// measurement comparison
// ********************

TCanvas* c_calicompthick[maxmeas];
TH2D* h_calicompthick[maxmeas];
TLegend* l_calicompthick[maxmeas];

TCanvas* c_calicompmaterial[maxmeas];
TH2D* h_calicompmaterial[maxmeas];
TLegend* l_calicompmaterial[maxmeas];


TCanvas* c_blockcompthick[maxmeas];
TH1D* h_blockcompthick[maxmeas];
TLegend* l_blockcompthick[maxmeas];

TCanvas* c_blockcompmaterial[maxmeas];
TH1D* h_blockcompmaterial[maxmeas];
TLegend* l_blockcompmaterial[maxmeas];
TH2D* h_blockcompmaterial_2D[maxmeas];

TGraphErrors* g_blockcompmaterial[maxmeas];
TCanvas* c_blockcompmaterial_g;
TH2D* h_blockcompmaterial_g;
TLegend* l_blockcompmaterial_g;

TGraphErrors* g_blockcompmaterial2[maxmeas];
TCanvas* c_blockcompmaterial_g2;
TH2D* h_blockcompmaterial_g2;
TLegend* l_blockcompmaterial_g2;

TGraphErrors* g_gradcompmaterial[maxmeas];
TCanvas* c_gradcompmaterial;
TH2D* h_gradcompmaterial;
TLegend* l_gradcompmaterial;

TH2D* h_calicomp[sensors];


// ********************
// this function prepares the output file and all the roots
// ********************

void prepareroot()
{

	if (debug<5)
	{
		cout << " " << endl;
		cout << "********************" << endl;
		cout << "Booking ROOTs!" << endl;
		cout << "********************" << endl;
		cout << " " << endl;
	}

	// first all roots for each individual measurement

	// loop over the max number of measurements
	for (int ii=0;ii<maxmeas;ii++)
	{

		// a char for naming things
		char tempchar[100];

		// calibration canvas
		sprintf(tempchar, "c_cali%i", ii);
		c_cali[ii] = new TCanvas(tempchar,"Calibrations",600,400);
		c_cali[ii]->SetGrid();
		c_cali[ii]->SetFillColor(0);
		c_cali[ii]->SetBorderMode(0);
		c_cali[ii]->SetBorderSize(2);
		c_cali[ii]->SetFrameBorderMode(0);
		c_cali[ii]->SetFrameBorderMode(0);

		// graph for the average calibration in a measurement run
		avg_calibrationgraph[ii] = new TGraphErrors();

		// the histogram for calibrations
		sprintf(tempchar, "h_cali%i", ii);
		h_cali[ii]= new TH2D(tempchar,"Calibrations", 9 , -0.5, 9.5, 100 , -10 , 10);
		h_cali[ii]->SetXTitle("Sensor");
		h_cali[ii]->SetYTitle("Calibration [#circC]");
		h_cali[ii]->SetStats(0000);

		// the legend for calibrations
		l_cali[ii] = new TLegend(0.7,0.1,0.9,0.3);
		l_cali[ii]->SetBorderSize(1);
		l_cali[ii]->SetFillColor(0);
		l_cali[ii]->SetFillStyle(1);


		// the canvas for temperatures
		sprintf(tempchar, "c_temps%i", ii);
		c_temps[ii]= new TCanvas(tempchar,"Temperatures",600,400);
		c_temps[ii]->SetGrid();
		c_temps[ii]->SetFillColor(0);
		c_temps[ii]->SetBorderMode(0);
		c_temps[ii]->SetBorderSize(2);
		c_temps[ii]->SetFrameBorderMode(0);
		c_temps[ii]->SetFrameBorderMode(0);

		// the histogram for temperatures
		sprintf(tempchar, "h_temps%i", ii);
		h_temps[ii] = new TH2D(tempchar,"Temperatures", 100 , 0, 120000, 100 , -10 , 50);
		h_temps[ii]->SetXTitle("Time [s]");
		h_temps[ii]->SetYTitle("Temperature [#circC]");
		h_temps[ii]->SetStats(0000);

		// the legend for temperatures
		l_temps[ii] = new TLegend(0.7,0.1,0.9,0.3);
		l_temps[ii]->SetBorderSize(1);
		l_temps[ii]->SetFillColor(0);
		l_temps[ii]->SetFillStyle(1);


		// the canvas for delta temperatures
		sprintf(tempchar, "c_deltatemps%i", ii);
		c_deltatemps[ii] = new TCanvas(tempchar,"Delta Temperatures",600,400);
		c_deltatemps[ii]->SetGrid();
		c_deltatemps[ii]->SetFillColor(0);
		c_deltatemps[ii]->SetBorderMode(0);
		c_deltatemps[ii]->SetBorderSize(2);
		c_deltatemps[ii]->SetFrameBorderMode(0);
		c_deltatemps[ii]->SetFrameBorderMode(0);

		// the histogram for delta temperatures
		sprintf(tempchar, "h_deltatemps%i", ii);
		h_deltatemps[ii] = new TH2D(tempchar,"Delta Temperatures", 100 , 0, 120000, 100 , -10 , 10);
		h_deltatemps[ii]->SetXTitle("Time [s]");
		h_deltatemps[ii]->SetYTitle("Delta Temperature [#circC]");
		h_deltatemps[ii]->SetStats(0000);

		// the legend for delta temperatures
		l_deltatemps[ii] = new TLegend(0.7,0.1,0.9,0.3);
		l_deltatemps[ii]->SetBorderSize(1);
		l_deltatemps[ii]->SetFillColor(0);
		l_deltatemps[ii]->SetFillStyle(1);


		// the canvas for gradients
		sprintf(tempchar, "c_gradtemps%i", ii);
		c_gradtemps[ii] = new TCanvas(tempchar,"Temperature Gradients",600,400);
		c_gradtemps[ii]->SetGrid();
		c_gradtemps[ii]->SetFillColor(0);
		c_gradtemps[ii]->SetBorderMode(0);
		c_gradtemps[ii]->SetBorderSize(2);
		c_gradtemps[ii]->SetFrameBorderMode(0);
		c_gradtemps[ii]->SetFrameBorderMode(0);

		// the histogram for gradients
		sprintf(tempchar, "h_gradtemps%i", ii);
		h_gradtemps[ii] = new TH2D(tempchar,"Temperature Gradients", 100 , 0, 80, 100 , -10 , 50);
		h_gradtemps[ii]->SetXTitle("Sensor Position [mm]");
		h_gradtemps[ii]->SetYTitle("Temperature [#circC]");
		h_gradtemps[ii]->SetStats(0000);

		// the legend for gradients
		l_gradtemps[ii] = new TLegend(0.7,0.1,0.9,0.3);
		l_gradtemps[ii]->SetBorderSize(1);
		l_gradtemps[ii]->SetFillColor(0);
		l_gradtemps[ii]->SetFillStyle(1);

		// the histogram for the difference between blocks
		sprintf(tempchar, "h_blockdifference%i", ii);
		h_blockdifference[ii] = new TH1D(tempchar,"Block Temperature Differences", 1000, -5, 10);
		h_blockdifference[ii]->SetXTitle("Temperature Difference [#circC]");
		h_blockdifference[ii]->SetYTitle("Entries");
		h_blockdifference[ii]->SetStats(1111);


		// the graphs for the individual calibrations
		for (int i=0;i<maxcalibs;i++)
		{
			calibrationgraph[ii][i] = new TGraphErrors();
			caliposition[ii][i] = new TLine();
		}


		// the graphs for the temperatures and delta temperatures of each individual sensor
		for (int i=0;i<sensors;i++)
		{
			tempgraph[ii][i] = new TGraphErrors();
			deltatempgraph[ii][i] = new TGraphErrors();
		}


		// the graphs for each gradient measurement point
		for (int i=0;i<maxstable;i++)
		{
			gradgraph[ii][i] = new TGraphErrors();
			sprintf(tempchar, "fitlow%i%i", ii,i);
			gradfit1[ii][i] = new TF1(tempchar, "pol1", 0.0, 40.1);
			sprintf(tempchar, "fithigh%i%i", ii,i);
			gradfit2[ii][i] = new TF1(tempchar, "pol1", 39.9, 80.0);
			gradposition[ii][i] = new TLine();
		}

	} // done measurement loop

	if (debug<2)
	{
		cout << "Done booking in measurement loop!" << endl;
	}

	// now plots for comparisons between measurements

	// for each thickness...
	for (int i=0;i<thicknesscount;i++)
	{

		// a char for naming things
		char tempchar[100];

		// calibrations
		sprintf(tempchar, "Calibrations at thickness %i", thicknesslist.at(i));
		c_calicompthick[i] = new TCanvas(tempchar,"Calibrations",600,400);
		c_calicompthick[i]->SetGrid();
		c_calicompthick[i]->SetFillColor(0);
		c_calicompthick[i]->SetBorderMode(0);
		c_calicompthick[i]->SetBorderSize(2);
		c_calicompthick[i]->SetFrameBorderMode(0);
		c_calicompthick[i]->SetFrameBorderMode(0);

		h_calicompthick[i]= new TH2D(tempchar,"Calibrations", 9 , -0.5, 9.5, 100 , -10 , 10);
		h_calicompthick[i]->SetXTitle("Sensor");
		h_calicompthick[i]->SetYTitle("Calibration [#circC]");
		h_calicompthick[i]->SetStats(0000);

		l_calicompthick[i] = new TLegend(0.7,0.1,0.9,0.3);
		l_calicompthick[i]->SetBorderSize(1);
		l_calicompthick[i]->SetFillColor(0);
		l_calicompthick[i]->SetFillStyle(1);


		// block temperature differences
		sprintf(tempchar, "Block temperature difference at thickness %i", thicknesslist.at(i));
		c_blockcompthick[i] = new TCanvas(tempchar,"Block temperature difference",600,400);
		c_blockcompthick[i]->SetGrid();
		c_blockcompthick[i]->SetFillColor(0);
		c_blockcompthick[i]->SetBorderMode(0);
		c_blockcompthick[i]->SetBorderSize(2);
		c_blockcompthick[i]->SetFrameBorderMode(0);
		c_blockcompthick[i]->SetFrameBorderMode(0);

		h_blockcompthick[i]= new TH1D(tempchar,"Block temperature difference", 1000, -5, 10);
		h_blockcompthick[i]->SetXTitle("Temperature [#circC]");
		h_blockcompthick[i]->SetYTitle("Entries");
		h_blockcompthick[i]->SetStats(0000);

		l_blockcompthick[i] = new TLegend(0.7,0.1,0.9,0.3);
		l_blockcompthick[i]->SetBorderSize(1);
		l_blockcompthick[i]->SetFillColor(0);
		l_blockcompthick[i]->SetFillStyle(1);

	}

	if (debug<2)
	{
		cout << "Done thickness loop!" << endl;
	}

	// for each material...
	for (int i=0;i<materialcount;i++)
	{

		// a char for naming things
		char tempchar[100];

		// calibrations
		sprintf(tempchar, "Calibrations of material %s", materiallist.at(i).c_str());
		c_calicompmaterial[i] = new TCanvas(tempchar,"Calibrations",600,400);
		c_calicompmaterial[i]->SetGrid();
		c_calicompmaterial[i]->SetFillColor(0);
		c_calicompmaterial[i]->SetBorderMode(0);
		c_calicompmaterial[i]->SetBorderSize(2);
		c_calicompmaterial[i]->SetFrameBorderMode(0);
		c_calicompmaterial[i]->SetFrameBorderMode(0);

		h_calicompmaterial[i]= new TH2D(tempchar,"Calibrations", 9 , -0.5, 9.5, 100 , -10 , 10);
		h_calicompmaterial[i]->SetXTitle("Sensor");
		h_calicompmaterial[i]->SetYTitle("Calibration [#circC]");
		h_calicompmaterial[i]->SetStats(0000);

		l_calicompmaterial[i] = new TLegend(0.7,0.1,0.9,0.3);
		l_calicompmaterial[i]->SetBorderSize(1);
		l_calicompmaterial[i]->SetFillColor(0);
		l_calicompmaterial[i]->SetFillStyle(1);


		// block temperature differences
		sprintf(tempchar, "Block temperature difference for material %s", materiallist.at(i).c_str());
		c_blockcompmaterial[i] = new TCanvas(tempchar,tempchar,600,400);
		c_blockcompmaterial[i]->SetGrid();
		c_blockcompmaterial[i]->SetFillColor(0);
		c_blockcompmaterial[i]->SetBorderMode(0);
		c_blockcompmaterial[i]->SetBorderSize(2);
		c_blockcompmaterial[i]->SetFrameBorderMode(0);
		c_blockcompmaterial[i]->SetFrameBorderMode(0);

		h_blockcompmaterial[i]= new TH1D(tempchar,tempchar, 1000, -5, 10);
		h_blockcompmaterial[i]->SetXTitle("Temperature [#circC]");
		h_blockcompmaterial[i]->SetYTitle("Entries");
		h_blockcompmaterial[i]->SetStats(1111);

		l_blockcompmaterial[i] = new TLegend(0.7,0.1,0.9,0.3);
		l_blockcompmaterial[i]->SetBorderSize(1);
		l_blockcompmaterial[i]->SetFillColor(0);
		l_blockcompmaterial[i]->SetFillStyle(1);

		sprintf(tempchar, "Block temperature difference for material %s vs temperature", materiallist.at(i).c_str());
		h_blockcompmaterial_2D[i]= new TH2D(tempchar,tempchar, 100, 0, 30, 1000, -5, 10);
		h_blockcompmaterial_2D[i]->SetXTitle("Measurement Temperature [#circC]");
		h_blockcompmaterial_2D[i]->SetYTitle("Temperature Difference [#circC]");
		h_blockcompmaterial_2D[i]->SetStats(1111);

		g_blockcompmaterial[i] = new TGraphErrors();
		g_blockcompmaterial2[i] = new TGraphErrors();

		g_gradcompmaterial[i] = new TGraphErrors();

	} // done material loop

	if (debug<2)
	{
		cout << "Done booking in material loop!" << endl;
	}

	char tempchar[100];

	// an overall comparison
	sprintf(tempchar, "Block temperature differences");
	c_blockcompmaterial_g = new TCanvas(tempchar,tempchar,600,400);
	c_blockcompmaterial_g->SetGrid();
	c_blockcompmaterial_g->SetFillColor(0);
	c_blockcompmaterial_g->SetBorderMode(0);
	c_blockcompmaterial_g->SetBorderSize(2);
	c_blockcompmaterial_g->SetFrameBorderMode(0);
	c_blockcompmaterial_g->SetFrameBorderMode(0);

	h_blockcompmaterial_g = new TH2D(tempchar,tempchar, 1000, 0, 30, 1000, -5, 10);
	h_blockcompmaterial_g->SetXTitle("Measurement Temperature [#circC]");
	h_blockcompmaterial_g->SetYTitle("Temperature Difference [#circC]");
	h_blockcompmaterial_g->SetStats(0000);

	l_blockcompmaterial_g = new TLegend(0.7,0.1,0.9,0.3);
	l_blockcompmaterial_g->SetBorderSize(1);
	l_blockcompmaterial_g->SetFillColor(0);
	l_blockcompmaterial_g->SetFillStyle(1);


	// comparison vs slope
	sprintf(tempchar, "Block temperature differences vs slope");
	c_blockcompmaterial_g2 = new TCanvas(tempchar,tempchar,600,400);
	c_blockcompmaterial_g2->SetGrid();
	c_blockcompmaterial_g2->SetFillColor(0);
	c_blockcompmaterial_g2->SetBorderMode(0);
	c_blockcompmaterial_g2->SetBorderSize(2);
	c_blockcompmaterial_g2->SetFrameBorderMode(0);
	c_blockcompmaterial_g2->SetFrameBorderMode(0);

	h_blockcompmaterial_g2 = new TH2D(tempchar,tempchar, 1000, -0.01, 0.05, 1000, -5, 10);
	h_blockcompmaterial_g2->SetXTitle("Gradient Slope [#circC/mm]");
	h_blockcompmaterial_g2->SetYTitle("Temperature Difference [#circC]");
	h_blockcompmaterial_g2->SetStats(0000);

	l_blockcompmaterial_g2 = new TLegend(0.7,0.1,0.9,0.3);
	l_blockcompmaterial_g2->SetBorderSize(1);
	l_blockcompmaterial_g2->SetFillColor(0);
	l_blockcompmaterial_g2->SetFillStyle(1);


	// an overall comparison of gradients
	sprintf(tempchar, "Gradient differences");
	c_gradcompmaterial = new TCanvas(tempchar,tempchar,600,400);
	c_gradcompmaterial->SetGrid();
	c_gradcompmaterial->SetFillColor(0);
	c_gradcompmaterial->SetBorderMode(0);
	c_gradcompmaterial->SetBorderSize(2);
	c_gradcompmaterial->SetFrameBorderMode(0);
	c_gradcompmaterial->SetFrameBorderMode(0);

	h_gradcompmaterial = new TH2D(tempchar,tempchar, 1000, 0, 30, 1000, -0.01, 0.05);
	h_gradcompmaterial->SetXTitle("Heat - Workpoint [#circC]");
	h_gradcompmaterial->SetYTitle("Gradient Slope [#circC/mm]");
	h_gradcompmaterial->SetStats(0000);

	l_gradcompmaterial = new TLegend(0.7,0.1,0.9,0.3);
	l_gradcompmaterial->SetBorderSize(1);
	l_gradcompmaterial->SetFillColor(0);
	l_gradcompmaterial->SetFillStyle(1);

	// an overall comparison for calibrations
	for (int i=0;i<sensors;i++)
	{
		sprintf(tempchar, "Overall Calibration of Sensor %i", i);
		h_calicomp[i] = new TH2D(tempchar,tempchar, 30, 0, 30, 150, -5, 10);
		h_calicomp[i]->SetXTitle("Measurement Temperature [#circC]");
		h_calicomp[i]->SetYTitle("Calibration [#circC]");
		h_calicomp[i]->SetStats(1111);
	}

	if (debug<2)
	{
		cout << "Done overall booking!" << endl;
	}
}


// ********************
// this function opens the individual root file
// ********************

void openfile(std::string myfile)
{

	// the stream - check if file exists
	ifstream fileRead;

	// open
	fileRead.open(myfile.c_str());

	if ( !fileRead.is_open() )
	{
		cout << "Error opening root file " << myfile << " !" << endl;
		exit ( EXIT_FAILURE );
	}

	// close, let root do this
	fileRead.close();
  
	// the file to open
	TFile *f1;
	f1 = TFile::Open(myfile.c_str());

	// enter the input file
	f1->cd();

	// go into the tree
	mytuple = (TTree*)f1->Get("thermoDAQ");

	// how many entries in this tuple?
	tupleentrycount = mytuple->GetEntries();
	if (debug<5)
	{
		cout << " " << endl;
		cout << "********************" << endl;
		cout << "Looping "<< tupleentrycount << " entries in the input file " << myfile.c_str() << " ..." << endl;
		cout << "********************" << endl;
		cout << " " << endl;
	}

	// ********************
	// connect branch and variable
	// ********************

	// the time
	mytuple->SetBranchAddress("uTime", &uTime);

	// the temperatures of the sensors
	for (int row = 0; row < 10; ++row)
	{
		mytuple->SetBranchAddress(Form("temperature%d", row), &temperature[row]);
	}

	// the current
	mytuple->SetBranchAddress("current1", &current1);

	// the working temperature
	mytuple->SetBranchAddress("workingTemperature", &workingTemperature);

}


// ********************
// read the runlist - this opens the runlist and gets the tuple names and other infos for plotting
// ********************

void readrunlist(std::string astring)
{
	// the stream
	ifstream fileRead;

	// open
	fileRead.open(astring.c_str());

	if ( !fileRead.is_open() )
	{
		cout << "Error opening runlist file " << astring << " !" << endl;
		exit ( EXIT_FAILURE );
	}

	// the line to read
	std::string line;

	// a counter
	int linecount = 0;

	if (debug<5)
	{
		cout << " " << endl;
		cout << "********************" << endl;
		cout << "Reading runlist!" << endl;
		cout << "********************" << endl;
		cout << " " << endl;
	}

	int filecounter = 0;

	// loop over the lines in the file
	while (std::getline(fileRead, line))
	{

		// comment lines start with #
		string startpart = line.substr(0,1);
		if (startpart != "#")
		{

			// dummy string
			std::string fail;

			// the input string of this line
			std::string input;
			std::istringstream iss(line);
			input = iss.str();

			// the delimiter of the comments (space)
			std::string delimiter = ",";

			// the position we found of the delimiter
			size_t pos = 0;

			// the iterator for counting elements
			std::vector<int>::iterator it;

			// string iterator
			std::vector<string>::iterator its;

			// first the filename
			pos = input.find(delimiter);
			fail = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());

			// push the name into the vector
			filelist.push_back(fail);

			if (debug<4)
			{
				cout << "Found file no: " << filecounter << " : " << fail << " !" << endl;
			}

			// nameing starts at 0
			filecounter++;

			// the sensor sorting
			pos = input.find(delimiter);
			fail = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());

			// push the name into the vector
			sensorsort.push_back(fail.c_str());

			if (debug<4)
			{
				cout << "Found sensor sorting: " << fail << " !" << endl;
			}
			
			// the sensor sorting
			pos = input.find(delimiter);
			fail = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());

			// push the name into the vector
			brokenlist.push_back(fail.c_str());

			if (debug<4)
			{
				cout << "Found broken sensor: " << fail << " !" << endl;
			}

			// the material thickness
			pos = input.find(delimiter);
			fail = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());

			// push the name into the vector
			thickness.push_back(atoi(fail.c_str()));

			if (debug<4)
			{
				cout << "Found thickness: " << fail << " !" << endl;
			}

			// sort thicknesses
			it = find (thicknesslist.begin(), thicknesslist.end(), atoi(fail.c_str()));
			if (it != thicknesslist.end())
			{
				if (debug <1)
				{
					cout << "Thickness " << fail << " already in list!" << endl;
				}
			} else {
				thicknesscount++;
				thicknesslist.push_back(atoi(fail.c_str()));
				if (debug <2)
				{
					cout << "Thickness " << fail << " not in list! Total thicknesses now " << thicknesscount << endl;
				}
			}

			// the material type
			pos = input.find(delimiter);
			fail = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());

			// push the name into the vector
			material.push_back(fail);

			if (debug<4)
			{
				cout << "Found material: " << fail << " !" << endl;
			}

			// sort materials
			its = find (materiallist.begin(), materiallist.end(), fail);
			if (its != materiallist.end())
			{
				if (debug <1)
				{
					cout << "Material " << fail << " already in list!" << endl;
				}
			} else {
				materialcount++;
				materiallist.push_back(fail);
				if (debug <2)
				{
					cout << "Material " << fail << " not in list! Total materials now " << materialcount << endl;
				}
			}

			// comments
			pos = input.find(delimiter);
			fail = input.substr(0, pos);
			input.erase(0, pos + delimiter.length());

			// push the name into the vector
			comments.push_back(fail);

			if (debug<4)
			{
				cout << "Found comment: " << fail << " !" << endl;
				cout << " " << endl;
			}

			linecount++;

		} // done no comment line

	} // done while reading

	fileRead.close();

	if (linecount>maxmeas)
	{
		cout << "Warning: found " << linecount << " measurements in the runlist. Increase maxmeas from " << maxmeas << " to " << linecount << " to analyse all runs!" << endl;
	}

}


// ********************
// a function to evaluate if the change in sensor temperatures is below the given limit -> are we in thermal equilibrium?
// ********************

bool stablesystem(double mydelta)
{
	bool systemstablility = true;
	for (int i=0;i<sensors;i++)
	{
		deltaT[i] = temperature[i] - temperaturebefore[i];
		temperaturebefore[i] = temperature[i];
		// skip non-connected sensors //FIXME
		if (i != 0 && i!= 9)
		{
			if (!(fabs(deltaT[i]) <= mydelta))
			{
				systemstablility = false;
				if (debug<1)
				{
					cout << "System not in thermal equilibrium: deltaT [ " << i << " ] = " << deltaT[i] << endl;
				}
			}
		}
	}
	return(systemstablility);
}


// ********************
// a function to return the time
// ********************

double nicetime()
{

	// get the time of this event
	TDatime thetime;
	thetime = TDatime(uTime);

	// problem with the time: e.g. 17:59:59 is when called as thetime.GetTime() 175.959
	// the next time is 180000 as in 18:00:00 O'clock! Jump of 4000
	// get the corrected, continuous time
	if(usedpoints == 1)
	{
		//this sets j to the next hour "thetime.GetTime()" will hit
		time0 = thetime.GetTime();
		temptime = (thetime.GetTime()-(thetime.GetTime() % 10000))/10000 + 1;
	}

	if(temptime == 24 && thetime.GetTime() < 1000)
	{
		  // every full day we jump 240.000 minus the 4000 for a full hour!
		  time0= time0 - 236000;
		  temptime = 1;
	}

	if(thetime.GetTime() >= (temptime*10000))
	{
		  // every full hour we jump 4000 
		  time0 = time0 +4000;
		  temptime++;
	}

	double sometime = thetime.GetTime() - time0;
	// now we have a continuous time!

	return sometime;

}


// ********************
// a function to sort the sensors according to their position
// ********************

void sensorsorting(int run)
{
	float temptemperature[sensors] = {0.0};
	for (int j=0;j<sensors;j++)
	{
		temptemperature[j] = temperature[j];
	}
	for (size_t i=0;i<sensorsort.at(run).size();i++)
	{
		stringstream astream;
		astream << sensorsort.at(run).at(i);
		int temp = 0;
		astream >> temp;
		temperature[i] = temptemperature[temp];
		if (debug<1)
		{
			cout << "Sorting sensors. Position " << i << " now has sensor " << sensorsort.at(run).at(i) << " with temperature " << temperature[i] << " C" << endl;
		}
	}
}


// ********************
// a function to return the new position of a (broken) sorted sensor
// ********************

int brokensorting(int run, int origsensor)
{
	int newsensor = -1;
	for (size_t i=0;i<sensorsort.at(run).size();i++)
	{
		stringstream astream;
		astream << sensorsort.at(run).at(i);
		int temp = 0;
		astream >> temp;
		if (temp == origsensor)
		{
			newsensor = i;
		}
	}
	return newsensor;
}


// ********************
// the main function
// ********************

int main(int argc, char** argv)
{

	// user inputs which runlist to read
	std::string astring = "fail";
	stringstream astream;
	if (argc>1)
	{
		astream << argv[1];
		astring = astream.str();
	} else {
		cout << "You did not specify a runlist! Please input a runlist now:" << endl;
		getline (cin, astring);
	}

	// read the runlist into the vectors
	readrunlist(astring);

	// then prepare the roots
	prepareroot();

	// loop over each measurement in the runlist
	for (unsigned int ii=0;ii<filelist.size();ii++)
	{

		// open the file
		openfile(filelist.at(ii));

		// let's go!

		// prepare output
		outputFile->cd();
		char namechar[100];
		sprintf(namechar, "Measurement %i", ii);
		TDirectory* thisdirectory = outputFile->mkdir(namechar);
		thisdirectory->cd();


		// mode selection, 1 = calibration, 3 = calibration and analysis
		if (mode == 1 || mode == 3)
		{

			time1 = 0.0;
			temptime = 0.0;
			usedpoints = 0;

			// draw
			c_cali[ii]->cd();
			h_cali[ii]->Draw("");

			if (debug<5)
			{
				cout << " " << endl;
				cout << "********************" << endl;
				cout << "Running calibration!" << endl;
				cout << "********************" << endl;
				cout << " " << endl;
			}

			// loop over all tuple entries
			for(int i = 0; i < tupleentrycount;i++)
			{
				// only save every precision-th event
				if (i % precision == 0)
				{

					// to make sure there is a gap in between the calibrations, count points between
					inbetween++;
		
					// get the event
					mytuple->GetEntry(i);

					// count points for graphs
					usedpoints++;

					// set the time
					time1 = nicetime();

					// sort the sensors
					sensorsorting(ii);

					// now we analyse all calbration parts, there are 'calibs' of these.
					if(calibs[ii] < maxcalibs)
					{

						// there has to be a change of the bath temperature, then we start to look for a calibration point
						if (workingTemperature != workingTemperaturebefore && inbetween > 50)
						{
							lookForCal = true;
							if (debug<3)
							{
								cout << "Searching for calibration at time: " << time1 << endl;
								cout << "Current working temperature is: " << workingTemperature << endl;
								cout << "Measurement points since last search: " << inbetween << endl;
								cout << " " << endl;
							}
							inbetween=0;
						}

						// make sure there is a gap between the calibrations
						if(inbetween > 51)
						{
							workingTemperaturebefore = workingTemperature;
						}

						// the heater has to be off
						if (lookForCal && current1 == 0)
						{
							insidecool = true;
							if (debug<3)
							{
								cout << "Performing search!" << endl;
								cout << "Measurement time is: " << time1 << endl;
								cout << "Tuple point is: " << i << endl;
								cout << " " << endl;
							}

							insidecool = stablesystem(deltacali);

							// if alle sensors are good, then we add all temperatures to the calitemp
							if (insidecool)
							{
								if (debug<3)
								{
									cout << "All deltaTs are good!" << endl;
								}

								for (int j = 0; j < sensors ; j++)
								{
									calitemp[ii][j][calibs[ii]] = calitemp[ii][j][calibs[ii]] + temperature[j];
									if (debug<4)
									{
										cout << "Calibration point " << calitemp[ii][j][calibs[ii]] << " at working temperature " << workingTemperature << endl;
									}
									// shift the point for plotting
									calibrationgraph[ii][calibs[ii]]->SetPoint(j, j-(calibs[ii]/2 + 1)*0.1*pointsep,  calitemp[ii][j][calibs[ii]]-workingTemperature);
									calibrationgraph[ii][calibs[ii]]->SetPointError(j, 0,  (calitemp[ii][j][calibs[ii]]-workingTemperature)*errorpercentage);
								}
								calitime[ii][calibs[ii]] = time1;
								work_Temperature[ii][calibs[ii]]= workingTemperature;
								lookForCal = false;
								if (debug<4)
								{
									cout << " " << endl;
									cout << "Done " << calibs[ii] + 1 << " calibrations!" << endl;
									cout << " " << endl;
								}
								calibs[ii]++;
								pointsep = pointsep*(-1);
							}
						}
					} // done calibs loop
				} // done precision'th loop 
			} // done tuple loop

			// so now we can average the calitemps 
			if (calibs[ii] > 0)
			{

				if (debug<4)
				{
					cout << "Averaging calibration points of " << calibs[ii] << " calibrations!" <<endl;
				}

				float avg_cali_error[maxmeas][sensors] = {0.0};
				for (int j = 0; j < sensors; j++)
				{
					// get the average
					for (int k=0;k<calibs[ii];k++)
					{
						avg_calitemp[ii][j] += calitemp[ii][j][k]-work_Temperature[ii][k];
					}
					avg_calitemp[ii][j] /= calibs[ii];
					avg_calibrationgraph[ii]->SetPoint(j, j, avg_calitemp[ii][j]);

					// get the standard deviation
					for (int k=0;k<calibs[ii];k++)
					{
						avg_cali_error[ii][j] += (avg_calitemp[ii][j] - (calitemp[ii][j][k]-work_Temperature[ii][k]))*(avg_calitemp[ii][j] - (calitemp[ii][j][k]-work_Temperature[ii][k]));
					}
					avg_cali_error[ii][j] /= calibs[ii];
					avg_cali_error[ii][j] = sqrt(avg_cali_error[ii][j]);
					avg_calibrationgraph[ii]->SetPointError(j, 0, avg_cali_error[ii][j]);
					if (debug<4)
					{
						cout << "Average calibration of sensor " << j << " is " << avg_calitemp[ii][j] << " +- " << avg_cali_error[ii][j] << " deg C." << endl;
					}
				}
				if (debug<4)
				{
					cout << " " << endl;
				}

			// if no calibrations were found, we have a bad measurement!
			} else {
				if (debug<5)
				{
					cout << "No possible calibration point found in run " << ii << "!" << endl;
					cout << "Aborting run!" << endl;
					cout << " " << endl;
				}
				continue;
			}
			
			
			// plot the output
			for (int j=0;j<calibs[ii];j++)
			{
				c_cali[ii]->cd();
				calibrationgraph[ii][j]->SetMarkerStyle(34);
				calibrationgraph[ii][j]->SetMarkerColor(j+1);
				calibrationgraph[ii][j]->SetMarkerSize(2);
				calibrationgraph[ii][j]->SetLineColor(j+1);
				calibrationgraph[ii][j]->SetLineWidth(2);
				calibrationgraph[ii][j]->SetLineStyle(1);
				calibrationgraph[ii][j]->Draw("P");
				char tempchar[100];
				sprintf(tempchar, "%.1f #circC", work_Temperature[ii][j]);
				l_cali[ii]->AddEntry(calibrationgraph[ii][j],tempchar,"lp");
				c_cali[ii]->Update();
			}
			c_cali[ii]->cd();
			avg_calibrationgraph[ii]->SetMarkerStyle(34);
			avg_calibrationgraph[ii]->SetMarkerColor(calibs[ii]+1);
			avg_calibrationgraph[ii]->SetMarkerSize(2);
			avg_calibrationgraph[ii]->SetLineColor(calibs[ii]+1);
			avg_calibrationgraph[ii]->SetLineWidth(2);
			avg_calibrationgraph[ii]->SetLineStyle(1);
			avg_calibrationgraph[ii]->Draw("P");
			l_cali[ii]->AddEntry(avg_calibrationgraph[ii],"Average","lp");
			l_cali[ii]->Draw();
			c_cali[ii]->Update();
			c_cali[ii]->Write();
			c_cali[ii]->Close();

		} // done mode selection

		// mode selection, 2 = analysis, 3 = calibration and analysis
		if (mode == 2 || mode == 3)
		{

			if (debug<5)
			{
				cout << " " << endl;
				cout << "********************" << endl;
				cout << "Running analysis!" << endl;
				cout << "********************" << endl;
				cout << " " << endl;
			}

			time1 = 0.0;
			temptime = 0.0;
			usedpoints = 0;

			for (int j=0;j<sensors;j++)
			{
				temperaturebefore[j] = 0.0;
			}
			insidecool = false;
			inbetween = 0;

			// loop over all tuple entries
			for(int i = 0; i < tupleentrycount;i++)
			{
				// only save every precision-th event
				if (i % precision == 0)
				{
		
					// get the event
					mytuple->GetEntry(i);

					// count points for graphs
					usedpoints++;

					// set the time
					time1 = nicetime();

					// sort the sensors
					sensorsorting(ii);

					// apply calibration
					// use average if no calibration for a working point is found
					bool applyaverage = true;

					// if there is a calibration for this specific working temperature, apply it
					for (int j = 0;j<maxcalibs;j++)
					{
						// if the working temperature is within 5% of the one used for calibration
						if ((work_Temperature[ii][j] >= (workingTemperature-workingTemperature*0.05)) && (work_Temperature[ii][j] <= (workingTemperature+workingTemperature*0.05)))
						{
							temperature[0] = temperature[0] - calitemp[ii][0][j] + work_Temperature[ii][j];
							temperature[1] = temperature[1] - calitemp[ii][1][j] + work_Temperature[ii][j];
							temperature[2] = temperature[2] - calitemp[ii][2][j] + work_Temperature[ii][j];
							temperature[3] = temperature[3] - calitemp[ii][3][j] + work_Temperature[ii][j];
							temperature[4] = temperature[4] - calitemp[ii][4][j] + work_Temperature[ii][j];

							temperature[5] = temperature[5] - calitemp[ii][5][j] + work_Temperature[ii][j];
							temperature[6] = temperature[6] - calitemp[ii][6][j] + work_Temperature[ii][j];
							temperature[7] = temperature[7] - calitemp[ii][7][j] + work_Temperature[ii][j];
							temperature[8] = temperature[8] - calitemp[ii][8][j] + work_Temperature[ii][j];
							temperature[9] = temperature[9] - calitemp[ii][9][j] + work_Temperature[ii][j];

							if (debug<1)
							{
								cout << "Found correct calibration at point " << j << " with " << work_Temperature[ii][j] << " deg C!" << endl;
							}

							applyaverage = false;
							break;
						}
					}

					if (applyaverage)
					{
						temperature[0] = temperature[0] - avg_calitemp[ii][0];
						temperature[1] = temperature[1] - avg_calitemp[ii][1];
						temperature[2] = temperature[2] - avg_calitemp[ii][2];
						temperature[3] = temperature[3] - avg_calitemp[ii][3];
						temperature[4] = temperature[4] - avg_calitemp[ii][4];

						temperature[5] = temperature[5] - avg_calitemp[ii][5];
						temperature[6] = temperature[6] - avg_calitemp[ii][6];
						temperature[7] = temperature[7] - avg_calitemp[ii][7];
						temperature[8] = temperature[8] - avg_calitemp[ii][8];
						temperature[9] = temperature[9] - avg_calitemp[ii][9];

						if (debug<1)
						{
							cout << "Did not find correct calibration, applying average!" << endl;
						}
					}

					// the temperature difference between different sensors
					for (int j = 0;j<sensors;j++)
					{
						if (j > 0)
						{
							deltaDT[j] = temperature[j] - temperature[j-1];
						} else {
							deltaDT[j] = temperature[j] - temperature[9];
						}
					}

					// print some output
					if (debug<0)
					{
						cout << "Data point: " << i << " , temperatures: " << temperature[0] << " " << temperature[1] << " " << temperature[2] << " " << temperature[3] << " " << temperature[4] << " " << temperature[5] << " " << temperature[6] << " " << temperature[7] << " " << temperature[8] << " " << temperature[9] <<" at time: " << time1 <<  endl;
						cout << "DeltaDT is: " << deltaDT[0] << " " << deltaDT[1] << " " << deltaDT[2] << " " << deltaDT[3] << " " << deltaDT[4] << " " << deltaDT[5] << " " << deltaDT[6] << " " << deltaDT[7] << " " << deltaDT[8] << " " << deltaDT[9] << endl;
					}

					// fill the time graphs
					for (int j=0;j<sensors;j++)
					{
						tempgraph[ii][j]->SetPoint(usedpoints, time1, temperature[j]);
						deltatempgraph[ii][j]->SetPoint(usedpoints, time1, deltaDT[j]);
					}

					// are we in thermal equilibrium?
					insidecool = stablesystem(deltagrad);
					if (insidecool)
					{
						inbetween++;
					}

					// require >2 stable points between actual points, also current on
					if (insidecool && (inbetween > 2) && current1 > 0.0)
					{
						for (int j = 0; j < sensors ; j++)
						{
							stabletemp[ii][j][stablepoints[ii]] = temperature[j];
							
						}
						if (debug<3)
						{
							cout << "Found stable point no. " << stablepoints[ii] << " at " << time1 << " s!" << endl;
						}

						// save working temperature, current and time
						stablework[ii][stablepoints[ii]] = workingTemperature;
						stablecurrent[ii][stablepoints[ii]] = current1;
						stabletime[ii][stablepoints[ii]] = time1;

						// increase the count
						stablepoints[ii]++;

						// reset the distance counter
						inbetween = 0;
					}

				} // done precision'th loop 
			} // done tuple loop

			if (debug<4)
			{
				cout << "Done tuple loop, plotting output!" << endl;
				cout << "Found " << stablepoints[ii] << " gradient points!" << endl;
				cout << " " << endl;
			}

			// plot the output
			int tempcounter = 0;

			// the temperatures
			c_temps[ii]->cd();
			h_temps[ii]->Draw("");
			for (int j=0;j<sensors;j++)
			{

				// skip 0 //FIXME
				if (j != 0)
				{
					c_temps[ii]->cd();
					tempgraph[ii][j]->SetMarkerStyle(34);
					tempgraph[ii][j]->SetMarkerColor(tempcounter+1);
					tempgraph[ii][j]->SetMarkerSize(2);
					tempgraph[ii][j]->SetLineColor(tempcounter+1);
					tempgraph[ii][j]->SetLineWidth(2);
					tempgraph[ii][j]->SetLineStyle(1);
					tempgraph[ii][j]->Draw("L");
					char tempchar[100];
					sprintf(tempchar, "Sensor %i", j);
					l_temps[ii]->AddEntry(tempgraph[ii][j],tempchar,"lp");
					c_temps[ii]->Update();
					tempcounter++;
				}
			}
			tempcounter = 0;

			// draw the lines of the calibration times
			for (int j=0;j<calibs[ii];j++)
			{
				c_temps[ii]->cd();
				caliposition[ii][j]->SetLineWidth(1);
				caliposition[ii][j]->SetLineColor(2);
				caliposition[ii][j]->SetLineStyle(6);
				caliposition[ii][j]->SetX1(calitime[ii][j]);
				caliposition[ii][j]->SetX2(calitime[ii][j]);
				caliposition[ii][j]->SetY1(-10);
				caliposition[ii][j]->SetY2(50);
				caliposition[ii][j]->Draw();
			}

			// draw the lines of the gradient times
			for (int j=0;j<stablepoints[ii];j++)
			{
				c_temps[ii]->cd();
				gradposition[ii][j]->SetLineWidth(1);
				gradposition[ii][j]->SetLineColor(1);
				gradposition[ii][j]->SetLineStyle(5);
				gradposition[ii][j]->SetX1(stabletime[ii][j]);
				gradposition[ii][j]->SetX2(stabletime[ii][j]);
				gradposition[ii][j]->SetY1(-10);
				gradposition[ii][j]->SetY2(50);
				gradposition[ii][j]->Draw();
			}

			c_temps[ii]->cd();
			l_temps[ii]->Draw();
			c_temps[ii]->Update();
			c_temps[ii]->Write();
			c_temps[ii]->Close();

			// the delta temperatures
			c_deltatemps[ii]->cd();
			h_deltatemps[ii]->Draw("");
			for (int j=0;j<sensors;j++)
			{
				// skip 0 and 4
				// also 1 and 5
				if (j != 0 && j != 4 && j != 1 && j != 5)
				{
					c_deltatemps[ii]->cd();
					deltatempgraph[ii][j]->SetMarkerStyle(34);
					deltatempgraph[ii][j]->SetMarkerColor(tempcounter+1);
					deltatempgraph[ii][j]->SetMarkerSize(2);
					deltatempgraph[ii][j]->SetLineColor(tempcounter+1);
					deltatempgraph[ii][j]->SetLineWidth(2);
					deltatempgraph[ii][j]->SetLineStyle(1);
					deltatempgraph[ii][j]->Draw("L");
					char tempchar[100];
					sprintf(tempchar, "Sensor %i - Sensor %i", j, j-1);
					l_deltatemps[ii]->AddEntry(deltatempgraph[ii][j],tempchar,"lp");
					c_deltatemps[ii]->Update();
					tempcounter++;
				}
			}
			tempcounter = 0;

			// draw the lines of the gradient times
			for (int j=0;j<stablepoints[ii];j++)
			{
				c_deltatemps[ii]->cd();
				gradposition[ii][j]->SetLineWidth(1);
				gradposition[ii][j]->SetLineColor(1);
				gradposition[ii][j]->SetLineStyle(5);
				gradposition[ii][j]->SetX1(stabletime[ii][j]);
				gradposition[ii][j]->SetX2(stabletime[ii][j]);
				gradposition[ii][j]->SetY1(-10);
				gradposition[ii][j]->SetY2(10);
				gradposition[ii][j]->Draw();
			}

			c_deltatemps[ii]->cd();
			l_deltatemps[ii]->Draw();
			c_deltatemps[ii]->Update();
			c_deltatemps[ii]->Write();
			c_deltatemps[ii]->Close();

			// the temperature gradients
			c_gradtemps[ii]->cd();
			h_gradtemps[ii]->Draw();
			const int n = 8;
			for (int j=0;j<stablepoints[ii];j++)
			{
				Double_t x[n] = {72,64,56,48,32,24,16,8};
				Double_t y[n] = {stabletemp[ii][1][j],stabletemp[ii][2][j],stabletemp[ii][3][j],stabletemp[ii][4][j],stabletemp[ii][5][j],stabletemp[ii][6][j],stabletemp[ii][7][j],stabletemp[ii][8][j]};
				for (int k=0;k<8;k++)
				{

					if (debug<2)
					{
						cout << "Adding point " << k << " of stable point " << j << " at " << x[k] << " mm, " << y[k] << " K!" << endl;
					}
					gradgraph[ii][j]->SetPoint(k,x[k]-((j/2 + 1)*0.5*pointsep),y[k]);
					gradgraph[ii][j]->SetPointError(k,1,y[k]*errorpercentage);

				}
				pointsep = pointsep * (-1);
				c_gradtemps[ii]->cd();
				gradgraph[ii][j]->SetMarkerStyle(34);
				gradgraph[ii][j]->SetMarkerColor(j+1);
				gradgraph[ii][j]->SetMarkerSize(2);
				gradgraph[ii][j]->SetLineColor(j+1);
				gradgraph[ii][j]->SetLineWidth(2);
				gradgraph[ii][j]->SetLineStyle(1);
				gradgraph[ii][j]->Draw("P");
				char tempchar[100];
				sprintf(tempchar, "Measurement point %i", j);
				l_gradtemps[ii]->AddEntry(gradgraph[ii][j],tempchar,"lp");
				c_gradtemps[ii]->Update();
				
				// check if there is a bad sensor and remove it from the gradient plot:
				for (int k=0;k<sensors;k++)
				{
					std::ostringstream ss;
					ss << k;
					
					std::size_t found = brokenlist.at(ii).find(ss.str());
					if (found!=std::string::npos)
					{
						if (debug<3)
						{
							cout << "Found broken sensor " << k << " removing point " << brokensorting(ii,k) << endl;
						}
						gradgraph[ii][j]->RemovePoint(brokensorting(ii,k)-1);
					}
				}

				// fit
				sprintf(tempchar, "fitlow%i%i", ii, j);
				gradgraph[ii][j]->Fit(tempchar, "RQ");
				gradfit1[ii][j]->Draw("l same");
				sprintf(tempchar, "fithigh%i%i", ii, j);
				gradgraph[ii][j]->Fit(tempchar, "RQ");
				gradfit2[ii][j]->Draw("l same");

				// calculate the temperature difference from the fit difference
				float lowtemp = gradfit1[ii][j]->Eval(40.0) + greasetemp/2.0;
				float hightemp = gradfit2[ii][j]->Eval(40.0) - greasetemp/2.0;
				tempdiff[ii][j] = hightemp - lowtemp;

				// define the measurement temperature of this as the average between top and bottom blocks
				tempdifftemp[ii][j] = (hightemp + lowtemp)/2.0;
				avg_tempdiff += tempdiff[ii][j];
				h_blockdifference[ii]->Fill(tempdiff[ii][j]);

				if (debug<4)
				{
					cout << "Point " << j << ":" << endl;
					cout << "Temperature difference between blocks is " << tempdiff[ii][j] << " K." << endl;
				}

				// the gradient between the aluminium blocks
				gradient_blocks += gradfit1[ii][j]->GetParameter(1);
				gradient_blocks += gradfit2[ii][j]->GetParameter(1);

				// calculate lambda of the blocks
				float lambda_al = resistor * stablecurrent[ii][j] * stablecurrent[ii][j] / (( (gradfit1[ii][j]->GetParameter(1) + gradfit2[ii][j]->GetParameter(1)) / 2.0*1000.0) * area );
				if (debug<4)
				{
					cout << "Alu Lambda is " << lambda_al << " W/(mK) at T = " << stablework[ii][j] << " °C." <<  endl;
				}
				
				if (debug<4)
				{
					cout << "Thermal resistance is " << tempdiff[ii][j]/(resistor * stablecurrent[ii][j] * stablecurrent[ii][j]) << " K/W." <<  endl;
					cout << " " << endl;
				}

			}
			c_gradtemps[ii]->cd();
			l_gradtemps[ii]->Draw();
			c_gradtemps[ii]->Update();

			// calculate average temperature difference
			if (stablepoints[ii] > 0)
			{
				avg_tempdiff /= stablepoints[ii];
				if (debug<4)
				{
					cout << " " << endl;
					cout << "Average temperature difference between blocks is " << avg_tempdiff << " K!" << endl;
				}

				// the average gradient in aluminium, 2* since there are 2 fits...
				gradient_blocks /= (2*stablepoints[ii]);
				if (debug<4)
				{
					cout << "Average temperature gradient between blocks is " << gradient_blocks*1000.0 << " K/m." << endl;
					cout << " " << endl;
				}
			} else {
				if (debug<5)
				{
					cout << "No stable gradient points found!" << endl;
					cout << " " << endl;
				}
			}

			c_gradtemps[ii]->Write();
			c_gradtemps[ii]->Close();
			h_blockdifference[ii]->Write();

		} // done mode selection

		if (mode == 0)
		{
		  
			if (debug<5)
			{
				cout << " " << endl;
				cout << "********************" << endl;
				cout << "Running testing!" << endl;
				cout << "********************" << endl;
				cout << " " << endl;
			}

			// loop over all tuple entries
			for(int ij = 0; ij < tupleentrycount;ij++)
			{
				// only save every precision-th event
				if (ij % precision == 0)
				{
		
					// count points for graphs
					usedpoints++;
		
					//make sure, there is a gap in between the calibrations
					inbetween++;
		
					// get the event
					mytuple->GetEntry(ij);
		
					
				} // done precision'th loop 
			} // done tuple loop
		} // done mode selection

		// catch wrong mode entry
		if ( mode < 0 || mode > 3)
		{
			cout << "Wrong value for mode! Allowed settings are:" << endl;
			cout << " 0 - Testing" << endl;
			cout << " 1 - Calibration" << endl;
			cout << " 2 - Analysis" << endl;
			cout << " 3 - Calibration and analysis" << endl;
			continue;
		}

		// we're done with a file, so some cleaning up
		time1 = 0.0;
		temptime = 0.0;
		usedpoints = 0;

	} // done measurement loop

	// now time to do some comparisons between measurement runs

	if (debug<5)
	{
		cout << " " << endl;
		cout << "********************" << endl;
		cout << "Preparing final analysis!" << endl;
		cout << "********************" << endl;
		cout << " " << endl;
	}

	// point counter for this graph
	int g_blockcompmaterialcount[maxmeas] = {0};
	int g_gradcompmaterialcount[maxmeas] = {0};

	// prepare output canvas
	outputFile->cd();
	c_blockcompmaterial_g->cd();
	h_blockcompmaterial_g->Draw();

	c_blockcompmaterial_g2->cd();
	h_blockcompmaterial_g2->Draw();


	c_gradcompmaterial->cd();
	h_gradcompmaterial->Draw();

	// loop the materials
	for (int l=0;l<materialcount;l++)
	{

		if (debug<1)
		{
			cout << "Looping all materials: " << materiallist.at(l) << endl;
		}

		// loop thicknesses
		for (int m=0;m<thicknesscount;m++)
		{

			if (debug<1)
			{
				cout << "Looping all thicknesses: " << thicknesslist.at(m) << endl;
			}

			// go through all measurements and look for the different parameters - "sorting"
			for (unsigned int ik=0;ik<filelist.size();ik++)
			{

				if (material.at(ik) == materiallist.at(l))
				{
					if (thickness.at(ik) == thicknesslist.at(m))
					{
						if (debug<4)
						{
							cout << "Looping all calibrations!" << endl;
							cout << " " << endl;
						}

						// go over the calibrations
						for (int j=0;j<calibs[ik];j++)
						{
							for (int k=0;k<sensors;k++)
							{
								h_calicomp[k]->Fill(work_Temperature[ik][j], calitemp[ik][k][j]-work_Temperature[ik][j]);
								// some more comparison?
								if (debug<2)
								{
									cout << "Run: " << ik << ", calibration: " << j << " , sensor: " << k << " , calibration temperature: " << calitemp[ik][k][j] << " , working temperature: "<< work_Temperature[ik][j] << endl;
								}
							}
						}

						if (debug<2)
						{
							cout << " " << endl;
						}

						// go over the stable points
						for (int j=0;j<stablepoints[ik];j++)
						{
							// fill histograms
							// the temperature difference
							h_blockcompmaterial[l]->Fill( tempdiff[ik][j]);

							// the temperature difference at the measurement temp
							h_blockcompmaterial_2D[l]->Fill(tempdifftemp[ik][j] , tempdiff[ik][j]);

							// add the points to the graph of this material
							g_blockcompmaterial[l]->SetPoint(g_blockcompmaterialcount[l],tempdifftemp[ik][j],tempdiff[ik][j]);

							// vs slope

							// maxmeas maxstable
							float slope1 = gradfit1[ik][j]->GetParameter(1);
							float slope2 = gradfit2[ik][j]->GetParameter(1);
							g_blockcompmaterial2[l]->SetPoint(g_blockcompmaterialcount[l],((slope1 + slope2)/2.0),tempdiff[ik][j]);
							g_blockcompmaterialcount[l]++;

							float lambda = resistor * stablecurrent[ik][j] * stablecurrent[ik][j] / area / ((slope1 + slope2)/2.0*1000.0);

							if (debug<2)
							{
								cout << "Run: " << ik << ", stable point: " << j << ", lambda: " << lambda << endl;
							}

							g_gradcompmaterial[l]->SetPoint(g_gradcompmaterialcount[l],tempdifftemp[ik][j],lambda);
							
							g_gradcompmaterialcount[l]++;

							if (debug<1)
							{
								cout << "Filling comparison output histos: Run " << ik << ", stable point " << j << endl;
							}
						}

						if (debug<2)
						{
							cout << " " << endl;
						}

					} // done if thickness
				} // done if material

			} // done filelist loop int i
		
		} // done thickness loop int m

		// the comparisons between materials
		outputFile->cd();
		char tempchar[100];
		sprintf(tempchar, "Comparison of Material %s", materiallist.at(l).c_str());
		TDirectory* materialdirectory = outputFile->mkdir(tempchar);
		materialdirectory->cd();
		h_blockcompmaterial[l]->Write();
		h_blockcompmaterial_2D[l]->Write();

		outputFile->cd();
		c_blockcompmaterial_g->cd();
		g_blockcompmaterial[l]->SetMarkerStyle(34);
		g_blockcompmaterial[l]->SetMarkerColor(l+1);
		g_blockcompmaterial[l]->SetMarkerSize(2);
		g_blockcompmaterial[l]->SetLineColor(l+1);
		g_blockcompmaterial[l]->SetLineWidth(2);
		g_blockcompmaterial[l]->SetLineStyle(1);
		g_blockcompmaterial[l]->Draw("P");
		sprintf(tempchar, "Material %s", materiallist.at(l).c_str());
		l_blockcompmaterial_g->AddEntry(g_blockcompmaterial[l],tempchar,"lp");
		c_blockcompmaterial_g->Update();

		outputFile->cd();
		c_blockcompmaterial_g2->cd();
		g_blockcompmaterial2[l]->SetMarkerStyle(34);
		g_blockcompmaterial2[l]->SetMarkerColor(l+1);
		g_blockcompmaterial2[l]->SetMarkerSize(2);
		g_blockcompmaterial2[l]->SetLineColor(l+1);
		g_blockcompmaterial2[l]->SetLineWidth(2);
		g_blockcompmaterial2[l]->SetLineStyle(1);
		g_blockcompmaterial2[l]->Draw("P");
		sprintf(tempchar, "Material %s", materiallist.at(l).c_str());
		l_blockcompmaterial_g2->AddEntry(g_blockcompmaterial2[l],tempchar,"lp");
		c_blockcompmaterial_g2->Update();

		outputFile->cd();
		c_gradcompmaterial->cd();
		g_gradcompmaterial[l]->SetMarkerStyle(34);
		g_gradcompmaterial[l]->SetMarkerColor(l+1);
		g_gradcompmaterial[l]->SetMarkerSize(2);
		g_gradcompmaterial[l]->SetLineColor(l+1);
		g_gradcompmaterial[l]->SetLineWidth(2);
		g_gradcompmaterial[l]->SetLineStyle(1);
		g_gradcompmaterial[l]->Draw("P");
		sprintf(tempchar, "Material %s", materiallist.at(l).c_str());
		l_gradcompmaterial->AddEntry(g_gradcompmaterial[l],tempchar,"lp");
		c_gradcompmaterial->Update();

	} // done material loop int l
	
	c_blockcompmaterial_g->cd();
	l_blockcompmaterial_g->Draw();
	c_blockcompmaterial_g->Update();
	c_blockcompmaterial_g->Write();
	c_blockcompmaterial_g->Close();

	c_blockcompmaterial_g2->cd();
	l_blockcompmaterial_g2->Draw();
	c_blockcompmaterial_g2->Update();
	c_blockcompmaterial_g2->Write();
	c_blockcompmaterial_g2->Close();

	c_gradcompmaterial->cd();
	l_gradcompmaterial->Draw();
	c_gradcompmaterial->Update();
	c_gradcompmaterial->Write();
	c_gradcompmaterial->Close();
	
	outputFile->cd();
	TDirectory* calidirectory = outputFile->mkdir("Calibration Comparison");
	calidirectory->cd();
//	char tempchar[100];
	for (int i=0;i<sensors;i++)
	{
		h_calicomp[i]->Write();
	}
	outputFile->cd();

	// not sure if needed...
	return 0;

} // done main loop



