/**
 * @file NoiseGenerator.cxx
 * @brief Utility used to generate noise hits.
 * @author Omar Moreno, SLAC National Accelerator Laboratory
 */

#include "Tools/NoiseGenerator.h"
#include "Framework/Exception/Exception.h"
#include "Framework/EventProcessor.h"


#include <TH1F.h>
#include <TCanvas.h>
#include <vector>
#include <string>


namespace ldmx {

NoiseGenerator::NoiseGenerator(double noiseValue, bool gauss) {
  noise_ = noiseValue;
  useGaussianModel_ = gauss;
  poisson_dist_ =
      std::make_unique<boost::math::poisson_distribution<> >(noiseValue);
}

NoiseGenerator::~NoiseGenerator() {}

void NoiseGenerator::seedGenerator(uint64_t seed) {
  random_ = std::make_unique<TRandom3>(seed);
}

//void NoiseGenerator::onProcessStart(std::vector <std::vector<std::string>> histograms) {
//  getHistoDirectory();
//  for (const auto& histogramValues : histograms){
////    histograms_.create(histogramValues[0],histogramValues[1],std::stod (histogramValues[2]),std:stod(histogramValues[3]), std::stod(histogramValues[4]), histogramValues[5], std::stod(histogramValues[6]), std::stod(histogramValues[7]), std::stod(histogramValues[8]));
//    histograms_.create("RecoilX", "", 20, -0.5, 19.5, "RecoilX [mm]", 90, -450.0, 450.0);
//  }
//}
//    void NoiseGenerator::onProcessStart() {
//        farmework::getHistoDirectory();
//        histograms_.create("RecoilX", "", 20, -0.5, 19.5, "RecoilX [mm]", 90, -450.0, 450.0);
//    }

std::vector<double> NoiseGenerator::generateNoiseHits(int emptyChannels) {
  if (random_.get() == nullptr) {
    EXCEPTION_RAISE("RandomSeedException",
                    "Noise generator was not seeded before use");
  }
//   std::cout << "[ Noise Generator ]: Empty channels: "
//             << emptyChannels << std::endl;
//   std::cout << "[ Noise Generator ]: Normalized integration limit: "
//            << noiseThreshold_ << std::endl;

  double integral;
//  std::cout << "Use Gaussian bool value is: " << useGaussianModel_ << "; ";

//  std::cout << "[Poisson_dist, noise] is: " << complement(*poisson_dist_, noiseThreshold_ - 1) << "; ";
//  std::cout << "[noiseThreshold_ -1] is: " << noiseThreshold_ - 1  << "; ";
  double noise_Amp = 1;

  if (useGaussianModel_){
//    std::cout << "Noise  value is: " << noise_ << "; ";
//    std::cout << "pedestal value is: " << pedestal_ << std::endl;
    noise_ = noise_Amp* noise_;
    pedestal_ = 1 * pedestal_;
    noiseThreshold_ = 1 * noiseThreshold_;
    integral = ROOT::Math::normal_cdf_c(noiseThreshold_, noise_, pedestal_);}
  else //broken
    integral =
        boost::math::cdf(complement(*poisson_dist_, noiseThreshold_ - 1));
//   std::cout << "[ Noise Generator ]: Integral: "
//            << integral << std::endl;

  double noiseHitCount = random_->Binomial(emptyChannels, integral);
//   std::cout << "[ Noise Generator ]: # Noise hits: "
//            << noiseHitCount << std::endl;

    // Create a histogram
//    TH1F* hNoiseHitCounts = new TH1F("hNoiseHitCounts", "Noise Hit Counts", 10, 0, 10);
//
//    // Fill the histogram with the noise hit counts
//    hNoiseHitCounts->Fill(noiseHitCount);
//
//    // Create a canvas to draw the histogram
//    TCanvas* canvas = new TCanvas("canvas", "Histogram of Noise Hit Counts", 800, 600);
//    hNoiseHitCounts->Draw();
//
//    // Save the histogram to a file
//    canvas->SaveAs("/home/ananda/ldmx-results/noiseHitCounts"+ std::to_string(noise_Amp) +".png");
//
//    // Clean up
//    delete hNoiseHitCounts;
//    delete canvas;

  std::vector<double> noiseHits;
  for (int hitIndex = 0; hitIndex < noiseHitCount; ++hitIndex) {
    double rand = random_->Uniform();
//     std::cout << "[ Noise Generator ]: Rand: "
//              << rand << std::endl;
    double draw = integral * rand;
//     std::cout << "[ Noise Generator ]: Draw: "
//              << draw << std::endl;

    double cumulativeProb = 1.0 - integral + draw;
//     std::cout << "[ Noise Generator ]: Cumulative probability: "
//              << cumulativeProb << std::endl;

    double valueAboveThreshold;
    if (useGaussianModel_)
      valueAboveThreshold =
          ROOT::Math::gaussian_quantile(cumulativeProb, noise_);
    else
      valueAboveThreshold =
          boost::math::quantile(*poisson_dist_, cumulativeProb);
//     std::cout << "[ Noise Generator ]: Noise value: "
//              << valueAboveThreshold << std::endl;

    noiseHits.push_back(valueAboveThreshold);
  }

  return noiseHits;
}

}  // namespace ldmx
