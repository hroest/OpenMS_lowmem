// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include "MzMLConsumer.h"
#include "MSDataReader.h"

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>

#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

using namespace OpenMS;
using namespace std;

template <typename MapType>
class OPENMS_DLLAPI GaussMzMLConsumer :
  public Internal::MzMLConsumer <MapType>
{

public:

  GaussMzMLConsumer(String filename, GaussFilter gf) :
    Internal::MzMLConsumer<MapType>(filename, ProgressLogger())
  {
    gf_ = gf;
  }

  void processSpectrum_(typename MapType::SpectrumType & s)
  {
    gf_.filter(s);
  }

  GaussFilter gf_;
};

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
    @page UTIL_LowMemNoiseFilterGaussian LowMemNoiseFilterGaussian

    @brief  Executes a Gaussian filter to reduce the noise in an MS experiment.
    It uses less memory than the regular NoiseFilterGaussian but is somewhat
    slower.

<CENTER>
    <table>
        <tr>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. predecessor tools </td>
      <td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ LowMemNoiseFilterGaussian \f$ \longrightarrow \f$</td>
            <td ALIGN = "center" BGCOLOR="#EBEBEB"> pot. successor tools </td>
        </tr>
        <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FileConverter </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerWavelet</td>
        </tr>
        <tr>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_Resampler </td>
      <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_PeakPickerHiRes</td>
    </tr>
    <tr>
            <td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_BaselineFilter</td>
        </tr>
    </table>
</CENTER>

    The Gaussian filter is a peak area preserving low-pass filter and is characterized by narrow bandwidths,
    sharp cutoffs, and low passband ripple.

    @note The Gaussian filter works for uniform as well as for non-uniform data.

    <B>The command line parameters of this tool are:</B>
  @verbinclude TOPP_NoiseFilterGaussian.cli
    <B>INI file documentation of this tool:</B>
    @htmlinclude TOPP_NoiseFilterGaussian.html

    <B>The algorithm parameters for the Gaussian filter are:</B>
@htmlinclude OpenMS_GaussFilter.parameters

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPNoiseFilterGaussianLowMem :
  public TOPPBase
{
public:
  TOPPNoiseFilterGaussianLowMem() :
    TOPPBase("LowMemNoiseFilterGaussian", "Removes noise from profile spectra by using Gaussian filter (on uniform as well as non-uniform data).", false)
  {
  }

  void registerOptionsAndFlags_()
  {
    registerInputFile_("in", "<file>", "", "input raw data file ");
    setValidFormats_("in", StringList::create("mzML"));
    registerOutputFile_("out", "<file>", "", "output raw data file ");
    setValidFormats_("out", StringList::create("mzML"));

    registerSubsection_("algorithm", "Algorithm parameters section");
  }

  Param getSubsectionDefaults_(const String & /*section*/) const
  {
    return GaussFilter().getDefaults();
  }

  ExitCodes main_(int, const char **)
  {
    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");

    ///////////////////////////////////
    // Load the file only to count the number of spectra/chromatograms
    ///////////////////////////////////
    Internal::MSDataCounterReserve<Peak1D> exp_cnt;
    {
      MzMLFile f;
      f.getOptions().addMSLevel(-1);
      f.load(in, exp_cnt);
    }
    ///////////////////////////////////
    // Load the file again to get the experimental settings (e.g.
    // InstrumentSettings etc.) to store in the target file.
    ///////////////////////////////////
    MSExperiment<Peak1D> settings_only;
    {
      MzMLFile f;
      f.getOptions().addMSLevel(-1);
      f.load(in, settings_only);
    }

#if 0
    //-------------------------------------------------------------
    // loading input
    //-------------------------------------------------------------
    MzMLFile mz_data_file;
    mz_data_file.setLogType(log_type_);
    MSExperiment<Peak1D> exp;
    mz_data_file.load(in, exp);
    if (exp.empty() && exp.getChromatograms().size() == 0)
    {
      LOG_WARN << "The given file does not contain any conventional peak data, but might"
                  " contain chromatograms. This tool currently cannot handle them, sorry.";
      return INCOMPATIBLE_INPUT_DATA;
    }
    //check for peak type (profile data required)
    if (!exp.empty() && PeakTypeEstimator().estimateType(exp[0].begin(), exp[0].end()) == SpectrumSettings::PEAKS)
    {
      writeLog_("Warning: OpenMS peak type estimation indicates that this is not profile data!");
    }

    //check if spectra are sorted
    for (Size i = 0; i < exp.size(); ++i)
    {
      if (!exp[i].isSorted())
      {
        writeLog_("Error: Not all spectra are sorted according to peak m/z positions. Use FileFilter to sort the input!");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }

    //check if chromatograms are sorted
    for (Size i = 0; i < exp.getChromatograms().size(); ++i)
    {
      if (!exp.getChromatogram(i).isSorted())
      {
        writeLog_("Error: Not all chromatograms are sorted according to peak m/z positions. Use FileFilter to sort the input!");
        return INCOMPATIBLE_INPUT_DATA;
      }
    }
#endif

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    Param filter_param = getParam_().copy("algorithm:", true);
    writeDebug_("Parameters passed to filter", filter_param, 3);

    ///////////////////////////////////
    // Create GaussFilter and hand it to the GaussMzMLConsumer
    ///////////////////////////////////
    GaussFilter gauss;
    gauss.setLogType(log_type_);
    gauss.setParameters(filter_param);
    boost::shared_ptr<GaussMzMLConsumer<MSExperiment<> > > gaussConsumer(
      new GaussMzMLConsumer<MSExperiment<> >(out, gauss));

    gaussConsumer->addDataProcessing(getProcessingInfo_(DataProcessing::SMOOTHING));

    ///////////////////////////////////
    // Set experimental settings in the consumer
    ///////////////////////////////////
    gaussConsumer->setExpectedSize(exp_cnt.spectraCounts, exp_cnt.chromatogramCounts);
    gaussConsumer->setExperimentalSettings(settings_only);

    ///////////////////////////////////
    // Create new MSDataReader and set our consumer
    ///////////////////////////////////
    Internal::MSDataReader<GaussMzMLConsumer<MSExperiment<> >,
        Peak1D, ChromatogramPeak> exp_reader;
    exp_reader.setConsumer(gaussConsumer);

    MzMLFile mz_data_file;
    mz_data_file.load(in, exp_reader);

    return EXECUTION_OK;
  }

};


int main(int argc, const char ** argv)
{
  TOPPNoiseFilterGaussianLowMem tool;
  return tool.main(argc, argv);
}


