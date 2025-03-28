/*######################################################################################################################
#
#  RothC C++ version
#
#  This C++ version was translated from the Python code by Caterina Toscano, Antonio Volta and Enrico Balugani 03/03/2025
#
#  The Rothamsted Carbon Model: RothC
#  Developed by David Jenkinson and Kevin Coleman
#
#  INPUTS:
#
#  clay:  clay content of the soil (units: %)
#  depth: depth of topsoil (units: cm)
#  IOM: inert organic matter (t C /ha)
#  nsteps: number of timesteps
#
#  year:    year
#  month:   month (1-12)
#  modern:   %modern
#  TMP:      Air temperature (C)
#  Rain:     Rainfall (mm)
#  Evap:     open pan evaporation (mm)
#  C_inp:    carbon input to the soil each month (units: t C /ha)
#  FYM:      Farmyard manure input to the soil each month (units: t C /ha)
#  PC:       Plant cover (0 = no cover, 1 = covered by a crop)
#  DPM/RPM:  Ratio of DPM to RPM for carbon additions to the soil (units: none)
#
#  OUTPUTS:
#
#  All pools are carbon and not organic matter
#
#  DPM:   Decomposable Plant Material (units: t C /ha)
#  RPM:   Resistant Plant Material    (units: t C /ha)
#  Bio:   Microbial Biomass           (units: t C /ha)
#  Hum:   Humified Organic Matter     (units: t C /ha)
#  IOM:   Inert Organic Matter        (units: t C /ha)
#  SOC:   Soil Organic Matter / Total organic Matter (units: t C / ha)
#
#  DPM_Rage:   radiocarbon age of DPM
#  RPM_Rage:   radiocarbon age of RPM
#  Bio_Rage:   radiocarbon age of Bio
#  HUM_Rage:   radiocarbon age of Hum
#  Total_Rage: radiocarbon age of SOC (/ TOC)
#
#  SWC:       soil moisture deficit (mm per soil depth)
#  RM_TMP:    rate modifying fator for temperature (0.0 - ~5.0)
#  RM_Moist:  rate modifying fator for moisture (0.0 - 1.0)
#  RM_PC:     rate modifying fator for plant retainment (0.6 or 1.0)

######################################################################################################################*/


#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <cmath>
//#include <math.h>


#include "commonConstants.h"
#include "rothCplusplus.h"

using namespace std;

void Crit3D_RothCplusplusMaps::initialize(const gis::Crit3DRasterGrid& DEM)
{
    decomposablePlantMaterial->initializeGrid(DEM);
    resistantPlantMaterial->initializeGrid(DEM);
    microbialBiomass->initializeGrid(DEM);
    humifiedOrganicMatter->initializeGrid(DEM);
    inertOrganicMatter->initializeGrid(DEM);
    soilOrganicMatter->initializeGrid(DEM);
}
// Calculates the plant retainment modifying factor (RMF_PC)
double RMF_PC(bool PC) {
    double RM_PC;
    if (!PC) {
        RM_PC = 1.0;
    } else {
        RM_PC = 0.6;
    }
    return RM_PC;
}

// Calculates the rate modifying factor for moisture (RMF_Moist)
double RMF_Moist(double RAIN, double PEVAP, double clay, double depth, bool PC, double &SWC) {
    const double RMFMax = 1.0;
    const double RMFMin = 0.2;

    //calc soil water functions properties
    double SMDMax = -(20 + 1.3 * clay - 0.01 * (clay * clay));
    double SMDMaxAdj = SMDMax * depth / 23.0;
    double SMD1bar = 0.444 * SMDMaxAdj;
    double SMDBare = 0.556 * SMDMaxAdj;

    double DF = RAIN - 0.75 * PEVAP;

    double minSWCDF = std::min(0.0, SWC + DF);
    double minSMDBareSWC = std::min(SMDBare, SWC);
    if (PC) {
        SWC = std::max(SMDMaxAdj, minSWCDF);
    } else {
        SWC = std::max(minSMDBareSWC, minSWCDF);
    }
    double RM_Moist;
    if (SWC > SMD1bar) {
        RM_Moist = 1.0;
    } else {
        RM_Moist = RMFMin + (RMFMax - RMFMin) * (SMDMaxAdj - SWC) / (SMDMaxAdj - SMD1bar);
    }
    return RM_Moist;
}

double RMF_Moist(double monthlyBIC, double clay, double depth, bool PC, double &SWC) {
    const double RMFMax = 1.0;
    const double RMFMin = 0.2;

    //calc soil water functions properties
    double SMDMax = -(20 + 1.3 * clay - 0.01 * (clay * clay));
    double SMDMaxAdj = SMDMax * depth / 23.0;
    double SMD1bar = 0.444 * SMDMaxAdj;
    double SMDBare = 0.556 * SMDMaxAdj;

    double DF = monthlyBIC;

    double minSWCDF = std::min(0.0, SWC + DF);
    double minSMDBareSWC = std::min(SMDBare, SWC);
    if (PC) {
        SWC = std::max(SMDMaxAdj, minSWCDF);
    } else {
        SWC = std::max(minSMDBareSWC, minSWCDF);
    }
    double RM_Moist;
    if (SWC > SMD1bar) {
        RM_Moist = 1.0;
    } else {
        RM_Moist = RMFMin + (RMFMax - RMFMin) * (SMDMaxAdj - SWC) / (SMDMaxAdj - SMD1bar);
    }
    return RM_Moist;
}

// Calculates the rate modifying factor for temperature (RMF_Tmp)
double RMF_Tmp(double TEMP) {
    double RM_TMP;
    if (TEMP < -5.0) {
        RM_TMP = 0.0;
    } else {
        RM_TMP = 47.91 / (std::exp(106.06 / (TEMP + 18.27)) + 1.0);
    }
    return RM_TMP;
}

void decomp(int timeFact, double &decomposablePlantMatter, double &resistantPlantMatter, double &microbialBiomass, double &humifiedOrganicMatter, double &IOM, double &SOC, double &decomposablePlantMatter_Rage, double &resistantPlantMatter_Rage, double &microbialBiomass_Rage, double &humifiedOrganicMatter_Rage, double &IOM_Rage, double &Total_Rage, double &modernC, double &modifyingRate, double &clay, double &C_Inp, double &FYM_Inp, double &decomposablePlantMatter_resistantPlantMatter)
{
    const double decomposablePlantMatter_k = 10.0;
    const double resistantPlantMatter_k = 0.3;
    const double microbialBiomass_k = 0.66;
    const double humifiedOrganicMatter_k = 0.02;

    //const double conr = 0.0001244876401867718; // equivalent to std::log(2.0)/5568.0;
    double tstep = 1.0/timeFact; //monthly 1/12 or daily 1/365
    double exc = std::exp(-CONR*tstep);

    //decomposition
    double decomposablePlantMatter1 = decomposablePlantMatter*std::exp(-modifyingRate*decomposablePlantMatter_k*tstep);
    double resistantPlantMatter1 = resistantPlantMatter*std::exp(-modifyingRate*resistantPlantMatter_k*tstep);
    double microbialBiomass1 = microbialBiomass*std::exp(-modifyingRate*microbialBiomass_k*tstep);
    double humifiedOrganicMatter1 = humifiedOrganicMatter*std::exp(-modifyingRate*humifiedOrganicMatter_k*tstep);

    double decomposablePlantMatterDelta = decomposablePlantMatter - decomposablePlantMatter1;
    double resistantPlantMatterDelta = resistantPlantMatter - resistantPlantMatter1;
    double microbialBiomassDelta = microbialBiomass - microbialBiomass1;
    double humifiedOrganicMatterDelta = humifiedOrganicMatter - humifiedOrganicMatter1;

    //calculating redistribution of carbon into each pool
    double x = 1.67*(1.85+1.60*std::exp(-0.0786*clay));
    double ratioFactor[3];
    ratioFactor[0] = x / (x + 1);
    ratioFactor[1] = 0.46 / (x + 1);
    ratioFactor[2] = 0.54 / (x + 1);
    //proportion C from each pool into CO2, microbialBiomass and humifiedOrganicMatter
    double decomposablePlantMatterToCo2 = decomposablePlantMatterDelta * ratioFactor[0];
    double decomposablePlantMatterToMicrobialBiomass = decomposablePlantMatterDelta * ratioFactor[1];
    double decomposablePlantMatterToHumifiedOrganicMatter = decomposablePlantMatterDelta * ratioFactor[2];

    double resistantPlantMatterToCo2 = resistantPlantMatterDelta * ratioFactor[0];
    double resistantPlantMatterToMicrobialBiomass = resistantPlantMatterDelta * ratioFactor[1];
    double resistantPlantMatterToHumifiedOrganicMatter = resistantPlantMatterDelta * ratioFactor[2];

    double microbialBiomassToCo2 = microbialBiomassDelta * ratioFactor[0];
    double microbialBiomassToMicrobialBiomass = microbialBiomassDelta * ratioFactor[1];
    double microbialBiomassToHumifiedOrganicMatter = microbialBiomassDelta * ratioFactor[2];

    double humifiedOrganicMatter_co2 = humifiedOrganicMatterDelta * ratioFactor[0];
    double humifiedOrganicMatter_microbialBiomass = humifiedOrganicMatterDelta * ratioFactor[1];
    double humifiedOrganicMatter_humifiedOrganicMatter = humifiedOrganicMatterDelta * ratioFactor[2];

    //update C pools
    decomposablePlantMatter = decomposablePlantMatter1;
    resistantPlantMatter = resistantPlantMatter1;
    microbialBiomass = microbialBiomass1 + decomposablePlantMatterToMicrobialBiomass + resistantPlantMatterToMicrobialBiomass + microbialBiomassToMicrobialBiomass + humifiedOrganicMatter_microbialBiomass;
    humifiedOrganicMatter = humifiedOrganicMatter1 + decomposablePlantMatterToHumifiedOrganicMatter + resistantPlantMatterToHumifiedOrganicMatter + microbialBiomassToHumifiedOrganicMatter + humifiedOrganicMatter_humifiedOrganicMatter;

    //split plant C to decomposablePlantMatter and resistantPlantMatter
    double PI_C_decomposablePlantMatter = decomposablePlantMatter_resistantPlantMatter / (decomposablePlantMatter_resistantPlantMatter + 1.0) * C_Inp;
    double PI_C_resistantPlantMatter = 1.0 / (decomposablePlantMatter_resistantPlantMatter + 1.0) * C_Inp;

    //split FYM C to decomposablePlantMatter, resistantPlantMatter and humifiedOrganicMatter
    double FYM_C_decomposablePlantMatter = 0.49*FYM_Inp;
    double FYM_C_resistantPlantMatter = 0.49*FYM_Inp;
    double FYM_C_humifiedOrganicMatter = 0.02*FYM_Inp;

    //add plant C and FYM_C to decomposablePlantMatter, resistantPlantMatter and humifiedOrganicMatter
    decomposablePlantMatter = decomposablePlantMatter + PI_C_decomposablePlantMatter + FYM_C_decomposablePlantMatter;
    resistantPlantMatter = resistantPlantMatter + PI_C_resistantPlantMatter + FYM_C_resistantPlantMatter;
    humifiedOrganicMatter = humifiedOrganicMatter + FYM_C_humifiedOrganicMatter;

    //calc new ract of each pool
    double decomposablePlantMatter_Ract = decomposablePlantMatter1 * std::exp(-CONR*decomposablePlantMatter_Rage);
    double resistantPlantMatter_Ract = resistantPlantMatter1 * std::exp(-CONR*resistantPlantMatter_Rage);

    double microbialBiomass_Ract = microbialBiomass1 * std::exp(-CONR*microbialBiomass_Rage);
    double decomposablePlantMatter_microbialBiomass_Ract = decomposablePlantMatterToMicrobialBiomass * std::exp(-CONR*decomposablePlantMatter_Rage);
    double resistantPlantMatter_microbialBiomass_Ract = resistantPlantMatterToMicrobialBiomass * std::exp(-CONR*resistantPlantMatter_Rage);
    double microbialBiomass_microbialBiomass_Ract = microbialBiomassToMicrobialBiomass * std::exp(-CONR*microbialBiomass_Rage);
    double humifiedOrganicMatter_microbialBiomass_Ract = humifiedOrganicMatter_microbialBiomass * std::exp(-CONR*humifiedOrganicMatter_Rage);

    double humifiedOrganicMatter_Ract = humifiedOrganicMatter1 *std::exp(-CONR*humifiedOrganicMatter_Rage);
    double decomposablePlantMatter_humifiedOrganicMatter_Ract = decomposablePlantMatterToHumifiedOrganicMatter * std::exp(-CONR*decomposablePlantMatter_Rage);
    double resistantPlantMatter_humifiedOrganicMatter_Ract = resistantPlantMatterToHumifiedOrganicMatter * std::exp(-CONR*resistantPlantMatter_Rage);
    double microbialBiomass_humifiedOrganicMatter_Ract = microbialBiomassToHumifiedOrganicMatter * std::exp(-CONR*microbialBiomass_Rage);
    double humifiedOrganicMatter_humifiedOrganicMatter_Ract = humifiedOrganicMatter_humifiedOrganicMatter * std::exp(-CONR*humifiedOrganicMatter_Rage);

    double IOM_Ract = IOM * std::exp(-CONR*IOM_Rage);

    //assign new C from plant and FYM the correct age
    double PI_decomposablePlantMatter_Ract = modernC * PI_C_decomposablePlantMatter;
    double PI_resistantPlantMatter_Ract = modernC * PI_C_resistantPlantMatter;

    double FYM_decomposablePlantMatter_Ract = modernC * FYM_C_decomposablePlantMatter;
    double FYM_resistantPlantMatter_Ract = modernC * FYM_C_resistantPlantMatter;
    double FYM_humifiedOrganicMatter_Ract = modernC * FYM_C_humifiedOrganicMatter;

    // update ract for each pool
    double decomposablePlantMatter_Ract_new = FYM_decomposablePlantMatter_Ract + PI_decomposablePlantMatter_Ract + decomposablePlantMatter_Ract*exc;
    double resistantPlantMatter_Ract_new = FYM_resistantPlantMatter_Ract + PI_resistantPlantMatter_Ract + resistantPlantMatter_Ract*exc;

    double microbialBiomass_Ract_new = (microbialBiomass_Ract + decomposablePlantMatter_microbialBiomass_Ract + resistantPlantMatter_microbialBiomass_Ract + microbialBiomass_microbialBiomass_Ract + humifiedOrganicMatter_microbialBiomass_Ract )*exc;

    double humifiedOrganicMatter_Ract_new = FYM_humifiedOrganicMatter_Ract + (humifiedOrganicMatter_Ract + decomposablePlantMatter_humifiedOrganicMatter_Ract + resistantPlantMatter_humifiedOrganicMatter_Ract + microbialBiomass_humifiedOrganicMatter_Ract + humifiedOrganicMatter_humifiedOrganicMatter_Ract)*exc;

    SOC = decomposablePlantMatter + resistantPlantMatter + microbialBiomass + humifiedOrganicMatter + IOM;
    double Total_Ract = decomposablePlantMatter_Ract_new + resistantPlantMatter_Ract_new + microbialBiomass_Ract_new + humifiedOrganicMatter_Ract_new + IOM_Ract;

    //calculate rage of each pool
    if (decomposablePlantMatter <= EPSILON)
        decomposablePlantMatter_Rage = 0;
    else
        decomposablePlantMatter_Rage = (std::log(decomposablePlantMatter/decomposablePlantMatter_Ract_new) ) / CONR;


    if(resistantPlantMatter <= EPSILON)
        resistantPlantMatter_Rage = 0;
    else
        resistantPlantMatter_Rage = (std::log(resistantPlantMatter/resistantPlantMatter_Ract_new) ) / CONR;

    if(microbialBiomass <= EPSILON)
        microbialBiomass_Rage = 0;
    else
        microbialBiomass_Rage = ( std::log(microbialBiomass/microbialBiomass_Ract_new) ) / CONR;


    if(humifiedOrganicMatter <= EPSILON)
        humifiedOrganicMatter_Rage = 0;
    else
        humifiedOrganicMatter_Rage = ( std::log(humifiedOrganicMatter/humifiedOrganicMatter_Ract_new) ) / CONR;


    if(SOC <= EPSILON)
        Total_Rage = 0;
    else
        Total_Rage = ( std::log(SOC/Total_Ract) ) / CONR;

    return;
}

// The Rothamsted Carbon Model: RothC
void RothC(int timeFact, double &DPM, double &RPM, double &BIO, double &HUM, double &IOM, double &SOC, double &DPM_Rage, double &RPM_Rage, double &BIO_Rage, double &HUM_Rage, double &IOM_Rage, double &Total_Rage, double &modernC, double &clay, double &depth, double &TEMP, double &RAIN, double &WATERLOSS,bool isET0, bool &PC, double &DPM_RPM, double C_Inp, double FYM_Inp, double &SWC) {
    // Calculate RMFs
    double RM_TMP = RMF_Tmp(TEMP);
    double RM_Moist;
    if (isET0)
    {
        double monthlyBIC = RAIN - WATERLOSS;
        RM_Moist = RMF_Moist(monthlyBIC, clay, depth, PC, SWC);
    }
    else
    {
        RM_Moist = RMF_Moist(RAIN, WATERLOSS, clay, depth, PC, SWC);
    }

    double RM_PC = RMF_PC(PC);

    // Combine RMF's into one.
    double modifyingRate = RM_TMP * RM_Moist * RM_PC;

    decomp(timeFact, DPM, RPM, BIO, HUM, IOM, SOC, DPM_Rage, RPM_Rage, BIO_Rage, HUM_Rage, IOM_Rage, Total_Rage, modernC, modifyingRate, clay, C_Inp, FYM_Inp, DPM_RPM);

    return;
}



std::vector<std::vector<double>> leggi_csv(const std::string& nome_file) {
    std::vector<std::vector<double>> dati;
    std::ifstream file(nome_file);

    if (!file.is_open()) {
        std::cerr << "Errore nell'apertura del file: " << nome_file << std::endl;
        return dati;
    }

    std::string linea;
    std::getline(file, linea);

    while (std::getline(file, linea)) {
        std::vector<double> riga;
        std::istringstream ss(linea);
        std::string token;

        for (int i = 0; i < 10; ++i) {
            std::getline(ss, token, ',');
            try {
                riga.push_back(std::stod(token));
            } catch (const std::invalid_argument& e) {
                std::cerr << "Errore nella conversione di un valore in double: " << token << std::endl;
            }
        }
        if (!riga.empty())
            dati.push_back(riga);
    }
    file.close();
    return dati;
}

void scrivi_csv(const std::string& nome_file, const std::vector<std::vector<double>>& dati) {
    std::ofstream file(nome_file);

    if (!file.is_open()) {
        std::cerr << "Errore nell'apertura del file: " << nome_file << std::endl;
        return;
    }

    file << "index,Year,Month,DPM_t_C_ha,RPM_t_C_ha,BIO_t_C_ha,HUM_t_C_ha,IOM_t_C_ha,SOC_t_C_ha,deltaC" << std::endl;

    for (const auto& riga : dati) {
        std::stringstream ss;
        for (double valore : riga) {
            ss << valore << ",";
        }
        ss.seekp(-1, std::ios::end);
        ss << std::endl;
        file << ss.str();
    }

    file.close();
}

/*
int main()
{
    //set initial pool values
    double DPM = 0;
    double RPM = 0;
    double BIO = 0;
    double HUM = 0;
    double SOC = 0;

    double DPM_Rage = 0.0;
    double RPM_Rage = 0.0;
    double BIO_Rage = 0.0;
    double HUM_Rage = 0.0;
    double IOM_Rage = 50000.0;

    //set initial soil water content (deficit)
    double SWC = 0;
    double TOC1 = 0;

    //TODO: read in RothC input data file
    double clay = 13.0;     //[%]
    double depth = 25.0;    //[cm]
    double IOM = 3.0041;    //[t C/ha]
    int nsteps = 840;       //[-]

    //std::vector<std::vector<double>> data = createDataMatrix();
    std::vector<std::vector<double>> data = leggi_csv("C:/Github/rothCStandAlone/data_input.csv");

    int k = -1;
    int j = -1;

    SOC = DPM + RPM + BIO + HUM + IOM;

    std::cout << j << "," << DPM << ","<< RPM << ","<< BIO << ","<< HUM << ","<< IOM << ","<< SOC << "\n";

    int timeFact = 12;
    double Total_Rage;

    double TEMP;
    double RAIN;
    double PEVAP;
    bool isET0 = false;
    bool PC;
    double DPM_RPM;
    double C_Inp;
    double FYM_Inp;
    double modernC;

    double test = 100;
    while (test > 0.000001)
    {
        k = k+1;
        j = j+1;

        if (k == timeFact) k = 0;

        TEMP = data[k][3];
        RAIN = data[k][4];
        PEVAP = data[k][5];
        PC = bool(data[k][8]);
        DPM_RPM = data[k][9];
        C_Inp = data[k][6];
        FYM_Inp = data[k][7];
        modernC = data[k][2]/100;

        Total_Rage = 0;

        RothC(timeFact, DPM, RPM, BIO, HUM, IOM, SOC, DPM_Rage, RPM_Rage, BIO_Rage, HUM_Rage, IOM_Rage, Total_Rage, modernC, clay, depth, TEMP, RAIN, PEVAP, isET0, PC, DPM_RPM, C_Inp, FYM_Inp, SWC);

        if (((k+1)%timeFact) == 0)
        {
            double TOC0 = TOC1;
            TOC1 = DPM + RPM + BIO + HUM;
            test = fabs(TOC1-TOC0);
        }
    }

    double totalDelta = (std::exp(-Total_Rage/8035.0) - 1) * 1000;

    std::cout << j << "," << DPM << "," << RPM << "," << BIO << "," << HUM << "," << SOC << "," << totalDelta << "\n";

    std::vector<std::vector<double>> yearList;
//    std::vector<std::vector<double>> yearList = {{double(1), double(j+1), DPM, RPM, BIO, HUM, IOM, SOC, totalDelta}};


    std::vector<std::vector<double>> monthList;
    int timeFactIndex;

    for (int i = timeFact; i < nsteps; i++)
    {
        TEMP = data[i][3];
        RAIN = data[i][4];
        PEVAP = data[i][5];
        PC = bool(data[i][8]);
        DPM_RPM = data[i][9];
        C_Inp = data[i][6];
        FYM_Inp = data[i][7];
        modernC = data[i][2]/100;

        RothC(timeFact, DPM, RPM, BIO, HUM, IOM, SOC, DPM_Rage, RPM_Rage, BIO_Rage, HUM_Rage, IOM_Rage, Total_Rage, modernC, clay, depth, TEMP, RAIN, PEVAP, isET0, PC, DPM_RPM, C_Inp, FYM_Inp, SWC);

        totalDelta = (std::exp(-Total_Rage/8035.0) - 1.0) * 1000;

        //std::cout << C_Inp << "," << FYM_Inp << "," << TEMP << "," << RAIN << "," << PEVAP << "," << SWC << ","
                  //<< PC << "," << DPM <<"," << RPM <<"," << BIO <<"," << HUM <<"," << IOM <<"," << SOC << "\n";

        monthList.push_back({double(i-timeFact), double(data[i][0]), double(data[i][1]), DPM, RPM, BIO, HUM, IOM, SOC, totalDelta});

        if (int(data[i][1]) == timeFact)
        {
            timeFactIndex = int(i/timeFact);
            yearList.push_back({double(timeFactIndex), data[i][0], data[i][1], DPM, RPM, BIO, HUM, IOM, SOC, totalDelta});
            //std::cout << i << "," << DPM << "," << RPM << "," << BIO << "," << HUM << "," << IOM << "," << SOC << "," << totalDelta << "\n";
        }
    }

    scrivi_csv("C:/Github/rothCStandAlone/CMonthResults.csv", monthList);
    scrivi_csv("C:/Github/rothCStandAlone/CYearResults.csv", yearList);

    return 0;


}

*/




