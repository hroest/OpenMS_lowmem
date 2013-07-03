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

#ifndef OPENMS_LOWMEM_CACHEDMZMLFILEHANDLER
#define OPENMS_LOWMEM_CACHEDMZMLFILEHANDLER

#include "MSDataIOInterface.h"

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include "OpenMS/ANALYSIS/OPENSWATH/CachedmzML.h"
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <boost/shared_ptr.hpp>

#include "MSDataReader.h"
#include "MzMLConsumer.h"

namespace OpenMS
{
  /**
      @brief File adapter for caching MzML files

      This class allows to cache an MzML file to a temporary space and will
      return an access pointer to the spectra/chromatograms.

      This class is intended to be used for minimize the memory footprint by
      not loading a whole MSExperiment into memory but rather loading it
      spectrum by spectrum and immediately caching it to disk. The advantage of
      this might be reduced memory usage at the expensive of disk usage and
      execution time.

      @ingroup FileIO
  */
  class OPENMS_DLLAPI CachedMzMLFileHandler 
  {

  public:
    /// Default constructor
    CachedMzMLFileHandler() { }

    /// Default destructor
    ~CachedMzMLFileHandler() {}

    /**
        @brief Caches an mzML file at a temporary location (e.g. "/tmp/")

        Will take an .mzML file at location "in" and cache it to 
        location tmp + tmp_fname. The caller is responsible that the temporary
        location is valid and can be written to. 

        String tmp = "/tmp/";
        String tmp_fname = "wf_tmpfile_0.mzML" ;
        OpenSwath::SpectrumAccessPtr spectra_ptr = 
          OpenMS::CachedMzMLFileHandler().cacheFile(filename, tmp, tmp_fname);
    */
    OpenSwath::SpectrumAccessPtr cacheLoadFile(String in, String tmp, String tmp_fname)
    {
      String cached_file = tmp + tmp_fname + ".cached";
      String meta_file = tmp + tmp_fname;
      cacheFile(in, tmp, tmp_fname);

      // Now re-load the cached file at tmp + tmp_fname and get a
      // SpectrumAccessPtr to return to caller.
      MSExperiment<Peak1D> exp;
      MzMLFile mzmlfile;
      mzmlfile.load(meta_file, exp);
      return SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
    }

    /**
        @brief Reads and input mzML file and caches it the locations given in cached_file and meta_file

        Assumes that the file at location "in" is a valid mzML file.
    */
    void cacheFile(String in, String cached_file, String meta_file)
    {

      // Load the file only to count the number of spectra/chromatograms
      Internal::MSDataCounterReserve<Peak1D> exp_cnt;
      {
        MzMLFile f;
        f.getOptions().addMSLevel(-1);
        f.load(in, exp_cnt);
      }

      // Open a cacher file with the known amounts of spectra/chromatograms at
      // the temporary location "cached_file"
      boost::shared_ptr<Internal::CachedMzMLConsumer> cacher(new Internal::CachedMzMLConsumer(cached_file));
      cacher->openFile(exp_cnt.spectraCounts, exp_cnt.chromatogramCounts);

      // Create a new file acceptor
      Internal::MSDataReader<Internal::CachedMzMLConsumer, Peak1D, ChromatogramPeak> exp_reader;
      exp_reader.setConsumer(cacher);
      {
        // Load the infile and write data directly to "cached_file"
        MzMLFile mz_data_file;
        mz_data_file.load(in, exp_reader);
      }

      //
      // Write metadata (e.g. ExperimentalSettings, SpectrumSettings etc.)
      //

      // Copy the metadata into a new structure
      MSExperiment<Peak1D> experiment_metadata;
      experiment_metadata = (ExperimentalSettings)exp_reader;
      experiment_metadata.setSpectra(exp_reader.getSpectraSettings());
      experiment_metadata.setChromatograms(exp_reader.getChromatogramSettings());

      // set dataprocessing on each spectrum/chromatogram
      DataProcessing dp;
      std::set<DataProcessing::ProcessingAction> actions;
      actions.insert(DataProcessing::FORMAT_CONVERSION);
      dp.setProcessingActions(actions);
      dp.setMetaValue("cached_data", "true");
      for (Size i=0; i<experiment_metadata.size(); ++i)
      {
        experiment_metadata[i].getDataProcessing().push_back(dp);
      }
      std::vector<MSChromatogram<ChromatogramPeak> > chromatograms = experiment_metadata.getChromatograms();
      for (Size i=0; i<chromatograms.size(); ++i)
      {
        chromatograms[i].getDataProcessing().push_back(dp);
      }
      experiment_metadata.setChromatograms(chromatograms);

      // write out metadata
      cacher->writeMetadata(experiment_metadata, meta_file);
    }

  };

} // NS OpenMS

#endif
