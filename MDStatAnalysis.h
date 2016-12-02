/*
 *  MDStatAnalysis.h
 *  
 *
 *  Created by Gregory Lehaut on 26/11/13.
 *  Copyright 2013 LPC Caen, CNRS/IN2P3. All rights reserved.
 *
 */


#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_eigen.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdlib.h>


class MDStatAnalysisTool{
protected:
	std::string infilename;
	
	size_t n_variable;
	
	gsl_vector * vevent;
	
	gsl_vector * vmean;
	gsl_vector * vsigma;
	gsl_matrix * m_covariance;
	gsl_matrix * m_correlation;
	gsl_vector * vmin;
	gsl_vector * vmax;
	
	gsl_matrix * m_eigenvector;
	gsl_vector * v_eigenvalue;
	gsl_vector * v_inertie;
	gsl_vector * v_cumul;
	
	
public:
	MDStatAnalysisTool();
	MDStatAnalysisTool( std::string filename, size_t dim );
	~MDStatAnalysisTool();
	
	void setInfilename( std::string filename );
	void setDimension( size_t );
	std::string getInfilename( ) const { return( infilename ); }
	
	void treatFile();
	void clear();
	void init();
	
	void printMeanInformation();
	void printSigmaInformation();
	void printCovarianceMatrix();
	void printCorrelationMatrix();
	
	void printEigenVector();
	void printEigenValue();
	
	void makeACP();
	
	void showEigenVector();
	
	void projectData();
};


