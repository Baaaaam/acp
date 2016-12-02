/*
 *  MDStatAnalysis.cc
 *  
 *
 *  Created by Gregory Lehaut on 26/11/13.
 *  Copyright 2013 LPC Caen, CNRS/IN2P3. All rights reserved.
 *
 *  g++ $(gsl-config --cflags) ./MDStatAnalysis.cc -o ./MDStatAnalysis $(gsl-config --libs) 
 */

#include "MDStatAnalysis.h"


int main( int argc_, char ** argv_ ){
	if (argc_ != 3 ) {
		std::clog << "MDStatAnalysis needs the dimension and the filename" << std::endl
		<< " $ MDStatAnalysis dimension infilename " << std::endl;
		exit(0);
	}
	size_t dim = atoi(argv_[1]);
	std::string infilename( argv_[2] );
	MDStatAnalysisTool analysis( infilename, dim );
	
	analysis.treatFile();
	std::clog << "********************************************************" << std::endl;
	analysis.printMeanInformation();
	std::clog << "********************************************************" << std::endl;
	
	analysis.printSigmaInformation();
	std::clog << "********************************************************" << std::endl;
	
	analysis.printCovarianceMatrix();
	std::clog << "********************************************************" << std::endl;
	
	analysis.printCorrelationMatrix();
	std::clog << "********************************************************" << std::endl;
	analysis.makeACP();
	std::clog << "********************************************************" << std::endl;
	analysis.printEigenValue();
	std::clog << "********************************************************" << std::endl;
	analysis.printEigenVector();
	std::clog << "********************************************************" << std::endl;
	analysis.showEigenVector();
	analysis.projectData();
	
	
	return( 1 );
}

/************************************************************************************/

MDStatAnalysisTool::MDStatAnalysisTool(){
	vmean = NULL;
	vsigma = NULL;
	v_inertie = NULL;
	v_cumul = NULL;
	vmin = NULL;
	vmax = NULL;
	m_covariance = NULL;
	m_correlation = NULL;
	v_eigenvalue = NULL;
	m_eigenvector = NULL;
	vevent = NULL;
}

MDStatAnalysisTool::MDStatAnalysisTool( std::string filename, size_t dim ){
	vmean = NULL;
	vsigma = NULL;
	v_inertie = NULL;
	v_cumul = NULL;
	m_covariance = NULL;
	m_correlation = NULL;
	vmin = NULL;
	vmax = NULL;
	vevent = NULL;
	v_eigenvalue = NULL;
	m_eigenvector = NULL;
	infilename = filename;
	setDimension( dim );
}

MDStatAnalysisTool::~MDStatAnalysisTool(){
	clear();
}

/* */
void MDStatAnalysisTool::setInfilename( std::string filename ){ infilename = filename; }
void MDStatAnalysisTool::setDimension( size_t dim ){ n_variable = dim; init(); }
void MDStatAnalysisTool::init(){
	clear();
	vmean = gsl_vector_calloc( n_variable );
	vmin = gsl_vector_calloc( n_variable );
	vmax = gsl_vector_calloc( n_variable );
	v_eigenvalue = gsl_vector_calloc( n_variable );
	vsigma = gsl_vector_calloc( n_variable );
	vevent = gsl_vector_calloc( n_variable );
	m_covariance = gsl_matrix_calloc( n_variable, n_variable );
	m_correlation = gsl_matrix_calloc( n_variable, n_variable );
	m_eigenvector = gsl_matrix_calloc( n_variable, n_variable );
	v_inertie = gsl_vector_calloc( n_variable );
	v_cumul = gsl_vector_calloc( n_variable );
}
void MDStatAnalysisTool::clear(){
	if (vmean != NULL) {
		gsl_vector_free( vmean );
		vmean = NULL;
	}
	if (vmin != NULL) {
		gsl_vector_free( vmin );
		vmin = NULL;
	}
	if (vmax != NULL) {
		gsl_vector_free( vmax );
		vmax = NULL;
	}
	if (v_eigenvalue != NULL) {
		gsl_vector_free( v_eigenvalue );
		v_eigenvalue = NULL;
	}
	if (v_inertie != NULL) {
		gsl_vector_free( v_inertie );
		v_inertie = NULL;
	}
	if (v_cumul != NULL) {
		gsl_vector_free( v_cumul );
		v_cumul = NULL;
	}
	if (vevent != NULL) {
		gsl_vector_free( vevent );
		vevent = NULL;
	}
	if (vsigma != NULL) {
		gsl_vector_free( vsigma );
		vsigma = NULL;
	}
	if (m_covariance != NULL) {
		gsl_matrix_free( m_covariance );
		m_covariance = NULL;
	}
	
	if (m_eigenvector != NULL) {
		gsl_matrix_free( m_eigenvector );
		m_eigenvector = NULL;
	}
	if (m_correlation != NULL) {
		gsl_matrix_free( m_correlation );
		m_correlation = NULL;
	}
}
void MDStatAnalysisTool::treatFile(){
	std::ifstream infile( infilename.c_str() );
	int nevent = 0;
	double tmpdouble;
	while (infile >> tmpdouble) {
		gsl_vector_set( vevent, 0, tmpdouble );
		for (size_t i = 1; i < n_variable; i++) {
			infile >> tmpdouble;
			gsl_vector_set( vevent, i, tmpdouble );
		}
		if (nevent == 0) {
			for (size_t i = 0; i < n_variable; i++) {
				gsl_vector_set( vmin, i, gsl_vector_get( vevent, i ) );
				gsl_vector_set( vmax, i, gsl_vector_get( vevent, i ) );
			}
		}
		else {
			for (size_t i = 0; i < n_variable; i++) {
				if (gsl_vector_get( vmax, i ) < gsl_vector_get( vevent, i )) {
					gsl_vector_set( vmax, i, gsl_vector_get( vevent, i ) );
				}
				
				if (gsl_vector_get( vmin, i ) > gsl_vector_get( vevent, i )) {
					gsl_vector_set( vmin, i, gsl_vector_get( vevent, i ) );
				}
			}
		}

		for (size_t i = 0; i < n_variable; i++) {
			gsl_vector_set( vmean, i, gsl_vector_get( vmean, i )+gsl_vector_get( vevent, i ) );
			gsl_vector_set( vsigma, i, gsl_vector_get( vsigma, i )+gsl_vector_get( vevent, i )*gsl_vector_get( vevent, i ) );
			for (size_t j = 0; j < n_variable; j++) {
				gsl_matrix_set( m_covariance, i, j, gsl_matrix_get( m_covariance, i, j )+gsl_vector_get( vevent, i )*gsl_vector_get( vevent, j ) );
			}
		}
		nevent++;
	};
	
	for (size_t i = 0; i < n_variable; i++) {
		gsl_vector_set( vmean, i, gsl_vector_get( vmean, i )/double( nevent ) );
		gsl_vector_set( vsigma, i, gsl_vector_get(vsigma, i ) /double( nevent ) - 
					   gsl_vector_get( vmean, i )* gsl_vector_get( vmean, i ) );
		gsl_vector_set( vsigma, i, sqrt( gsl_vector_get( vsigma, i ) ) );
	}
	for (size_t i = 0; i < n_variable; i++) {
		for (size_t j = 0; j < n_variable; j++) {
			gsl_matrix_set( m_covariance, i, j, gsl_matrix_get( m_covariance, i, j )/double(nevent) - gsl_vector_get(vmean, i ) * gsl_vector_get(vmean, j ) );
		}
	}
	for (size_t i = 0; i < n_variable; i++) {
		for (size_t j = 0; j < n_variable; j++) {
			gsl_matrix_set( m_correlation, i, j, gsl_matrix_get( m_covariance, i, j )/(gsl_vector_get(vsigma, i)*gsl_vector_get(vsigma,j)) );
		}
	}
	
	
	
	infile.close();
}

void MDStatAnalysisTool::printMeanInformation(){
	std::clog << "MDStatAnalysisTool Information, Mean value :" << std::endl;
	for (size_t i = 0; i < n_variable; i++) {
		std::clog << gsl_vector_get( vmean, i ) << " \t";
	}
	std::clog << std::endl;
}
void MDStatAnalysisTool::printSigmaInformation(){
	std::clog << "MDStatAnalysisTool Information, Sigma value :" << std::endl;
	for (size_t i = 0; i < n_variable; i++) {
		std::clog << gsl_vector_get( vsigma, i ) << " \t";
	}
	std::clog << std::endl;
}

void MDStatAnalysisTool::printCovarianceMatrix(){
	std::clog << "MDStatAnalysisTool Information,  covariance matrix :" << std::endl;
	for (size_t i = 0; i < n_variable; i++) {
		for (size_t j = 0; j < n_variable; j++) {
			std::clog << gsl_matrix_get( m_covariance, i, j ) << "   \t";
		}
		std::clog << std::endl;
	}
}


void MDStatAnalysisTool::printCorrelationMatrix(){
	std::clog << "MDStatAnalysisTool Information,  correlation matrix :" << std::endl;
	for (size_t i = 0; i < n_variable; i++) {
		for (size_t j = 0; j < n_variable; j++) {
			std::clog << gsl_matrix_get( m_correlation, i, j ) << "     \t";
		}
		std::clog << std::endl;
	}
}

void MDStatAnalysisTool::printEigenVector(){
	std::clog << "MDStatAnalysisTool Information,   Eigen Vector :" << std::endl;
	for (size_t i = 0; i < n_variable; i++) {
		double sum = 0.;
		for (size_t j = 0; j < n_variable; j++) {
			std::clog << gsl_matrix_get( m_eigenvector, i, j ) << "     \t";
			sum += gsl_matrix_get( m_eigenvector, i, j )*gsl_matrix_get( m_eigenvector, i, j );
		}
		std::clog << std::endl;
	}
}

void MDStatAnalysisTool::printEigenValue(){
	std::clog << "MDStatAnalysisTool Information, Eigenvalue :" << std::endl;
	double sum = 0.;
	std::clog << "Eigen value : \t";
	for (size_t i = 0; i < n_variable; i++) {
		sum += gsl_vector_get( v_eigenvalue, i );
		std::clog << gsl_vector_get( v_eigenvalue, i ) << "   \t";
	}
	std::clog << std::endl;
	
	std::clog << "normalize : \t";
	for (size_t i = 0; i < n_variable; i++) {
		std::clog << gsl_vector_get( v_eigenvalue, i )/sum << "   \t";
	}
	std::clog << std::endl;
	double sum2 =0;
	std::clog << "cumul : \t";
	for (size_t i = 0; i < n_variable; i++) {
		sum2 += gsl_vector_get( v_eigenvalue, i )/sum;
		std::clog << sum2 << "  \t";
	}
	std::clog << std::endl;
}

void MDStatAnalysisTool::makeACP(){
	gsl_eigen_symmv_workspace * w = gsl_eigen_symmv_alloc(n_variable);
	gsl_eigen_symmv( m_correlation, v_eigenvalue, m_eigenvector, w );
	gsl_eigen_symmv_sort( v_eigenvalue, m_eigenvector, GSL_EIGEN_SORT_VAL_DESC );
	gsl_eigen_symmv_free(w);
	double sum = 0.;
	for (size_t i = 0; i < n_variable; i++) {
		sum += gsl_vector_get( v_eigenvalue, i );
	}
	for (size_t i = 0; i < n_variable; i++) {
		gsl_vector_set( v_inertie, i, gsl_vector_get( v_eigenvalue, i )/sum );
	}
	gsl_vector_set( v_cumul, 0, gsl_vector_get( v_inertie, 0 ) );
	for (size_t i = 1; i < n_variable; i++) {
		gsl_vector_set( v_cumul, i, gsl_vector_get( v_inertie, i )+ gsl_vector_get( v_cumul, i-1 ) );
	}
	
}

void MDStatAnalysisTool::showEigenVector(){
	FILE * cmdgnu;
	
	cmdgnu = popen("gnuplot \n", "w");
	
		
	
	for (size_t i = 0; i < n_variable; i++) {
		fprintf(cmdgnu, "set style arrow %u head filled size screen 0.025,30,45 ls %u\n", i+1, i+1);
		fprintf(cmdgnu, "set arrow from 0,0 to %lf,%lf as %u\n",  gsl_matrix_get(m_eigenvector, i, 0),gsl_matrix_get(m_eigenvector, i, 1)  ,i+1);
		fprintf(cmdgnu, "set label 'x_%u' at 1.2,%lf right\n", i+1, 1.4-0.1*double(i));
		fprintf(cmdgnu, "set arrow from 1.3,%lf to 1.45,%lf as %u\n",  1.4-0.1*double(i),1.4-0.1*double(i) ,i+1);
		
	}
	fprintf(cmdgnu, "set title \"Projection of variables in the principal plan (%lf %%), file: %s\"\n",  gsl_vector_get(v_cumul,1)*100., infilename.c_str());
	fprintf(cmdgnu, "set xlabel \"new x1 (%lf %%)\"\n", gsl_vector_get(v_inertie,0)*100. );
	fprintf(cmdgnu, "set ylabel \"new x2 (%lf %%)\"\n", gsl_vector_get(v_inertie,1)*100. );
	fprintf(cmdgnu, "set grid\n");
	
	fprintf(cmdgnu, "set parametric \n" );
	fprintf(cmdgnu, "set dummy t\n");
	fprintf(cmdgnu, "set trange [-pi:pi]\n");
	fprintf(cmdgnu, "plot [-pi:pi][-1.5:1.5][-1.5:1.5] cos(t),sin(t) t\"\" ls 1, 0.5*cos(t),0.5*sin(t) t\"\" ls 1  \n");
	
	fflush(cmdgnu);
	fflush(cmdgnu);
	std::clog << "Enter a character : ...";
	char toto;
	std::cin >> toto;
	
	
	pclose(cmdgnu);
}

void MDStatAnalysisTool::projectData(){
	FILE * cmdgnu;
	
	cmdgnu = popen("gnuplot \n", "w");
	
		
	for (size_t i = 0; i < n_variable; i++) {
		fprintf(cmdgnu, "set style arrow %u head filled size screen 0.025,30,45 ls %u\n", i+1, i+1);
		fprintf(cmdgnu, "set arrow from 0,0 to %lf,%lf as %u\n",  gsl_matrix_get(m_eigenvector, i, 0),gsl_matrix_get(m_eigenvector, i, 1)  ,i+1);
		fprintf(cmdgnu, "set label 'x_%u' at 1.2,%lf right\n", i+1, 1.4-0.1*double(i));
		fprintf(cmdgnu, "set arrow from 1.3,%lf to 1.45,%lf as %u\n",  1.4-0.1*double(i),1.4-0.1*double(i) ,i+1);
		
	}
	fprintf(cmdgnu, "set title \"Projection of data in the principal plan (%lf %%), file: %s\"\n", gsl_vector_get(v_cumul,1)*100., infilename.c_str());
	fprintf(cmdgnu, "set xlabel \"new x1 (%lf %%)\"\n", gsl_vector_get(v_inertie,0)*100. );
	fprintf(cmdgnu, "set ylabel \"new x2 (%lf %%)\"\n", gsl_vector_get(v_inertie,1)*100. );
	fprintf(cmdgnu, "set grid\n");
	fprintf(cmdgnu, "set parametric \n" );
	fprintf(cmdgnu, "set dummy t\n");
	fprintf(cmdgnu, "set trange [-pi:pi]\n");
	fprintf(cmdgnu, "plot [-pi:pi][-1.5:1.5][-1.5:1.5] cos(t),sin(t) t\"\" ls 1, 0.5*cos(t),0.5*sin(t) t\"\" ls 1, \"-\" u($1):($2) t\"\" \n");
	std::ifstream infile( infilename.c_str() );
	int nevent = 0;
	double tmpdouble;
	gsl_vector * newX = gsl_vector_calloc(n_variable);
	while (infile >> tmpdouble) {
		gsl_vector_set( vevent, 0, tmpdouble );
		for (size_t i = 1; i < n_variable; i++) {
			infile >> tmpdouble;
			gsl_vector_set( vevent, i, tmpdouble );
		}
		for (size_t i = 0; i < n_variable; i++) {
			gsl_vector_set( newX, i, 0.0);
		}
		for (size_t i = 0; i < n_variable; i++) {
			for (size_t j = 0; j < n_variable; j++) {
				gsl_vector_set( newX, i, gsl_vector_get( newX, i )+gsl_matrix_get(m_eigenvector, j, i )*(gsl_vector_get(vevent, j )-(gsl_vector_get(vmax,j)+gsl_vector_get(vmin,j))*0.5)/(gsl_vector_get(vmax,j)-gsl_vector_get(vmin,j)) );
			}
		}
		fprintf(cmdgnu, "%lf %lf\n", gsl_vector_get(newX, 0 ), gsl_vector_get(newX, 1 ) );
	//	std::clog << gsl_vector_get(newX, 0 ) << "  " << gsl_vector_get(newX, 1 ) << std::endl;
		
	};
	fprintf(cmdgnu, "e\n");
	
	fflush(cmdgnu);
	fflush(cmdgnu);
	std::clog << "Enter a character : ...";
	char toto;
	std::cin >> toto;
	
	
	pclose(cmdgnu);

}