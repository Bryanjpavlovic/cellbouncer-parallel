// demux_vcf_io_species.cpp
// Species-mode I/O helpers for quant3_contam_ap.
// load_contam_prof lives in demux_vcf_io.cpp (added by demux_parallel session).
// This file adds only the species-level functions.

#include <string>
#include <algorithm>
#include <vector>
#include <map>
#include <set>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdlib>
#include <htswrapper/robin_hood/robin_hood.h>
#include "common.h"
#include "demux_vcf_io.h"

// =============================================================================
// load_species_prior (Appendix E)
// =============================================================================
void load_species_prior(const std::string& filename,
                        std::map<std::string, double>& species_prof,
                        std::map<std::string, double>& species_prof_conc){

    std::ifstream inf(filename);
    if (!inf.is_open()){
        fprintf(stderr, "ERROR: cannot open species_prior file: %s\n", filename.c_str());
        exit(1);
    }

    std::string line;
    while (std::getline(inf, line)){
        if (line.empty()) continue;
        std::istringstream splitter(line);
        std::string species_label;
        double proportion = 0.0;
        double alpha = -1.0;
        int fld = 0;
        std::string field;
        while (std::getline(splitter, field, '\t')){
            if (fld == 0) species_label = field;
            else if (fld == 1) proportion = std::atof(field.c_str());
            else if (fld == 2) alpha = std::atof(field.c_str());
            fld++;
        }
        if (species_label.empty()) continue;
        species_prof[species_label] = proportion;
        if (alpha >= 0) species_prof_conc[species_label] = alpha;
    }
}

// =============================================================================
// dump_species_prof
// =============================================================================
void dump_species_prof(FILE* outf,
                       const std::map<std::string, double>& species_prof,
                       const std::map<std::string, double>& species_prof_conc){
    for (const auto& sp : species_prof){
        if (species_prof_conc.count(sp.first) > 0){
            fprintf(outf, "%s\t%.10f\t%.6f\n",
                sp.first.c_str(), sp.second, species_prof_conc.at(sp.first));
        } else {
            fprintf(outf, "%s\t%.10f\n", sp.first.c_str(), sp.second);
        }
    }
}
