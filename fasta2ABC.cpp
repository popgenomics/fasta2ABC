#include <cstdlib>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <sstream>
#include <cmath>
#include <algorithm>

void getSequences( std::vector <std::string> & seq, std::vector <std::string> & nam, size_t & nSeq, const std::string fasta, const float minN);
size_t countN(const std::string & stringChain);
void testSequences(const size_t & nA, const size_t & nB, const size_t & minNInd);
void writeFiles(const std::vector <std::string> & seqA, const std::vector <std::string> & namA, const size_t & nA, const std::vector <std::string> & seqB, const std::vector <std::string> & namB, const size_t & nB, const std::string & geneName);
void positions_of_N(const std::vector <std::string> & seqA, const std::vector <std::string> & seqB, std::vector <size_t> & posN);
void positions_of_cleaned(const std::vector <std::string> & seqA, const std::vector <std::string> & seqB, std::vector <size_t> & pos_cleaned);
void checkCommandLine(const int argc);

int main(int argc, char* argv []){

	checkCommandLine(argc);

	const std::string fastaA(argv[1]); // name of fasta file for species A
	const std::string fastaB(argv[2]); // name of fasta file for species B 
	const float minN = std::stof(argv[3]); // maximum proportion of N allowed to keep an individual
	const size_t minNInd = std::stoi(argv[4]); //  minimum number of individuals to keep a species
	const std::string geneName = argv[5]; // geneName 
	
	size_t i(0);
	std::vector <std::string> seqA;	
	std::vector <std::string> namA;	
	size_t nA(0);
	
	std::vector <std::string> seqB;	
	std::vector <std::string> namB;	
	size_t nB(0);
	
	getSequences(seqA, namA, nA, fastaA, minN);
	getSequences(seqB, namB, nB, fastaB, minN);

	testSequences(nA, nB, minNInd); // stop the program if not enough sequences after cleaning in any species

	writeFiles(seqA, namA, nA, seqB, namB, nB, geneName);
	
	return(0);
}

void getSequences( std::vector <std::string> & seq, std::vector <std::string> & nam, size_t & nSeq, const std::string fasta, const float minN ){
	int i(-1);
	std::string line;
	std::ifstream fastaFile( fasta.c_str() );
	
	std::vector <std::string> seq_tmp;
	std::vector <std::string> nam_tmp;
	std::vector <size_t> nN;


	if( fastaFile ){
		while( std::getline( fastaFile, line) ){
			if( line[0] == '>' ){
				nN.push_back(0);
				nam_tmp.push_back( line.erase(0, 1) );
				seq_tmp.push_back( "" );
				++i;
			}else{
				seq_tmp[i].append(line);
			}
		}
	}else{
		std::cout << std::endl << "\tThe file " << fasta << " was not found" << std::endl << std::endl;
		exit(EXIT_FAILURE);
	}

	for( i=0; i<seq_tmp.size(); ++i){
		nN[i] = countN(seq_tmp[i]);
		if( nN[i] < minN*seq_tmp[0].size() ){
			seq.push_back(seq_tmp[i]);
			nam.push_back(nam_tmp[i]);
			++nSeq;
		}
	}
}

size_t countN(const std::string & stringChain){
	size_t i(0);
	size_t res(0);
	
	for(i=0; i<stringChain.size(); ++i){
		if( stringChain[i] == 'N' ){
			++res;
		}
	}
	return(res);
}

void testSequences(const size_t & nA, const size_t & nB, const size_t & minNInd){
	if( nA < minNInd){
		if( nB < minNInd ){
			std::cout << "\nNot enough individuals " << "(" << nA << " and " << nB << ") in both species after cleaning Ns\n" << std::endl;
			exit(EXIT_FAILURE);
		}else{
			std::cout << "\nNot enough individuals in species A" << " (" << nA << ") after cleaning Ns\n" << std::endl;
			exit(EXIT_FAILURE);
		}
	}else{
		if( nB < minNInd ){
			std::cout << "\nNot enough individuals in species B" << " (" << nB << ") after cleaning Ns\n" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
}


void writeFiles(const std::vector <std::string> & seqA, const std::vector <std::string> & namA, const size_t & nA, const std::vector <std::string> & seqB, const std::vector <std::string> & namB, const size_t & nB, const std::string & geneName){
	size_t i(0);
	size_t ind(0); // loop over individuals
	size_t pos(0); // loop over positions
	size_t pos_i(0); // loop over positions
	size_t test_pol(0); // ==0 if monomorphic; ==1 if polymorphic
	
	std::string ancestral(""); // ancestral allele = first found allele
	std::string allele(""); // currently readen allele
	std::vector <size_t> segSites; // positions of segsites

	std::string msStyle_tmp = "";	
	std::vector <std::string> msStyle; // msStyle[SNP][individual]
	
	std::vector <size_t> clean_positions; // contains positions for N
	positions_of_cleaned(seqA, seqB, clean_positions); // get positions for N
	
	const size_t nCleanedSites = clean_positions.size();
	for(pos=0; pos<nCleanedSites; ++pos){
		pos_i = clean_positions[pos];
		msStyle_tmp = "";
		test_pol = 0;	
			
		// loop over species A
		for(ind=0; ind<nA; ++ind){
			allele = seqA[ind][pos_i];
			
			if( ind==0 ){
				test_pol = 0;
				ancestral = allele;
			}
			if( allele[0] == ancestral[0] ){
				msStyle_tmp.append("0");
			}else{
				test_pol = 1;
				msStyle_tmp.append("1");
			}
		}
		
		// loop over species B
		for(ind=0; ind<nB; ++ind){
			allele = seqB[ind][pos_i];
			
			if( allele[0] == ancestral[0] ){
				msStyle_tmp.append("0");
			}else{
				test_pol = 1;
				msStyle_tmp.append("1");
			}
		}

		// if position is polymorphic
		if( test_pol == 1 ){
			segSites.push_back(pos_i);
			msStyle.push_back(msStyle_tmp);
		}
			
	} // end of loop over position

	const size_t nSegSites( segSites.size() );
	const size_t nSitesTot( clean_positions.size() );
	const size_t nSitesTot_withN( seqA[0].size() );

	// geneName_info.txt
	// locusName     L_noN     nSegSite     nsamA     nsamB
	std::ostringstream oss1;
	oss1 << geneName << "_info.txt";
	std::string infoFile_name = oss1.str();
	std::ofstream infoFile( infoFile_name.c_str(), std::ios::out );
	infoFile << "locusName\tL_noN\tnSegSite\tnsamA\tnsamB" << std::endl;
	infoFile << geneName << "\t" << nSitesTot << "\t" << nSegSites << "\t" << nA << "\t" << nB << std::endl;
	infoFile.close();
	
	// geneName.ms
	// locusName     L_noN     nSegSite     nsamA     nsamB
	std::ostringstream oss2;
	oss2 << geneName << ".ms";
	std::string msFile_name = oss2.str();
	std::ofstream msFile( msFile_name.c_str(), std::ios::out );
	msFile << "./msnsam tbs 20 -t tbs -r tbs tbs -I 2 tbs tbs 0 -m 1 2 tbs -m 2 1 tbs -n 1 tbs -n 2 tbs -ej tbs 1 2 -eN tbs tbs\n3579 27011 59243\n\n//\n";
	msFile << "segsites: " << nSegSites << std::endl;
	if( nSegSites > 0 ){
		msFile << "positions :";
		for(i=0; i<nSegSites; ++i){
			msFile << " " << (1.0*segSites[i]/nSitesTot_withN);
		}
		msFile << std::endl;
		
		for(ind=0; ind<(nA+nB); ++ind){
			for(pos=0; pos<nSegSites; ++pos){
				msFile << msStyle[pos][ind];
			}
			msFile << std::endl;
		}
	}
	msFile.close();

	// geneName.fasta
	std::ostringstream oss3;
	oss3 << geneName << "_interspe.fasta";
	std::string fastaFile_name = oss3.str();
	std::ofstream fastaFile( fastaFile_name.c_str(), std::ios::out );

	for(ind=0; ind<nA; ++ind){
		fastaFile << ">" << namA[ind] << std::endl;
		for(pos=0; pos<nCleanedSites; ++pos){
			fastaFile << seqA[ind][clean_positions[pos]];
		}
		fastaFile << std::endl;
	}
	
	for(ind=0; ind<nB; ++ind){
		fastaFile << ">" << namB[ind] << std::endl;
		for(pos=0; pos<nCleanedSites; ++pos){
			fastaFile << seqB[ind][clean_positions[pos]];
		}
		fastaFile << std::endl;
	}
}


void positions_of_N(const std::vector <std::string> & seqA, const std::vector <std::string> & seqB, std::vector <size_t> & posN){
	size_t pos_i(0);
	size_t ind_i(0);
	int test_pos(-1);
	for( ind_i = 0; ind_i < seqA.size(); ++ind_i){
		for( pos_i = 0; pos_i < seqA[ind_i].size(); ++pos_i ){
			if( seqA[ind_i][pos_i] == 'N' ){
				test_pos = -1;
				test_pos = std::count( posN.begin(), posN.end(), pos_i);
				if( test_pos == 0 ){
					posN.push_back(pos_i);
				}
			}
		}
	}
	for( ind_i = 0; ind_i < seqB.size(); ++ind_i){
		for( pos_i = 0; pos_i < seqB[ind_i].size(); ++pos_i ){
			if( seqB[ind_i][pos_i] == 'N' ){
				test_pos = -1;
				test_pos = std::count( posN.begin(), posN.end(), pos_i);
				if( test_pos == 0 ){
					posN.push_back(pos_i);
				}
			}
		}
	}
}


void positions_of_cleaned(const std::vector <std::string> & seqA, const std::vector <std::string> & seqB, std::vector <size_t> & pos_cleaned){
	size_t pos_i(0);
	size_t ind_i(0);
	int test_N(0);
	std::string alignement;

	for( pos_i=0; pos_i<seqA[0].size(); ++pos_i){
		alignement = "";
		for( ind_i=0; ind_i<seqA.size(); ++ind_i){
			alignement.push_back(seqA[ind_i][pos_i]);
		}
		
		for( ind_i=0; ind_i<seqB.size(); ++ind_i){
			alignement.push_back(seqB[ind_i][pos_i]);
		}
		
		test_N = 0;
		test_N = std::count( alignement.begin(), alignement.end(), 'N');
		if( test_N == 0 ){
			pos_cleaned.push_back(pos_i);
		}
	}
}


void checkCommandLine(const int argc){
	if( argc != 6 ){
		std::cout << std::endl << " getIntergenic produces intergenic sequences from a gff and a fasta file." << std::endl;
		std::cout << " in the current version: the original fasta file is only for one contig, but gff can contains informations for all contigs" << std::endl;
		std::cout << " 5 arguments are needed:" << std::endl;
		std::cout << "\tname of the fasta file for species A (string)." << std::endl;
		std::cout << "\tname of the fasta file for species B (string)." << std::endl;
		std::cout << "\tmaximum proportion of N allowed to keep an individual (float)." << std::endl;
		std::cout << "\tminimum number of individuals to keep a species (integer)." << std::endl;
		std::cout << "\tgene's name used for output files names (string)." << std::endl;
		std::cout << "\t\t./fasta2ABC locA_ama.fas locA_txn.fas 0.25 10 locusA" << std::endl << std::endl;
		std::cout << "\tcamille.roux.1983@gmail.com (30/06/2017)" << std::endl << std::endl;
		exit(EXIT_FAILURE);
	}
}

