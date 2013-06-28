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

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include "OpenMS/ANALYSIS/OPENSWATH/CachedmzML.h"

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <boost/shared_ptr.hpp>

#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>

#include <boost/function.hpp>

namespace OpenMS
{

  namespace InternalCaching
  {

    /*
     * Provides the same interface as MSExperiment in order to be passed into mzMLHandler (for example).
     *
     * Only the interface used by the mzMLHandler is implemented here and this
     * class should be used as base class. The child class then implements those
     * functions which it actually expects to be used.
    */
    template <typename PeakT = Peak1D, typename ChromatogramPeakT = ChromatogramPeak>
    class OPENMS_DLLAPI MSDataIOInterfaceBase :
      public ExperimentalSettings
    {

  public:

      /// @name Base type definitions
      //@{
      /// Peak type
      typedef PeakT PeakType;
      /// Chromatogram peak type
      typedef ChromatogramPeakT ChromatogramPeakType;
      /// Spectrum Type
      typedef MSSpectrum<PeakType> SpectrumType;
      //@}

      inline SpectrumType& operator[] (Size n)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return spectra_[n];
      }

      inline const SpectrumType& operator[] (Size n) const
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return spectra_[n];
      }

      inline Size size() const
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return spectra_.size(); 
      }

      inline void reserve(Size /* s */)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      inline bool empty() const
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return spectra_.empty(); 
      }
      //@}

      inline void reserveSpaceChromatograms(Size /* s */)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      inline void reserveSpaceSpectra(Size /* s */)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      /// Constructor
      MSDataIOInterfaceBase() :
        ExperimentalSettings()
      {
        // Please do not instantiate this
        // throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      /// Resets all internal values
      void reset()
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      /// adds a spectra to the list
      void addSpectrum(const MSSpectrum<PeakT> & spectrum)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      /// adds a chromatogram to the list
      void addChromatogram(const MSChromatogram<ChromatogramPeakType> & chromatogram)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

      /// returns the spectra list
      const std::vector<MSSpectrum<PeakT> > & getSpectra() const
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return spectra_;
      }

      /// returns the spectra list
      std::vector<MSSpectrum<PeakT> > & getSpectra() 
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return spectra_;
      }

      /// sets the chromatogram list
      void setChromatograms(const std::vector<MSChromatogram<ChromatogramPeakType> > & chromatograms)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        chromatograms_ = chromatograms;
      }

      /// returns the chromatogram list
      const std::vector<MSChromatogram<ChromatogramPeakType> > & getChromatograms() const
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return chromatograms_;
      }

      /// returns a single chromatogram 
      MSChromatogram<ChromatogramPeakType> & getChromatogram(Size id) 
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
        return chromatograms_[id];
      }

      /**
        @brief Clears all data and meta data

        @param clear_meta_data If @em true, all meta data is cleared in addition to the data.
      */
      void clear(bool clear_meta_data)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }

  protected:

      /// chromatograms
      std::vector<MSChromatogram<ChromatogramPeakType> > chromatograms_;
      /// spectra
      std::vector<SpectrumType> spectra_;
    };

    /*
     * Provides the same interface as MSExperiment in order to be passed into mzMLHandler (for example).
     *
     * Only the interface used by the mzMLHandler is implemented here and this
     * class should be used as base class. The child class then implements those
     * functions which it actually expects to be used.
     *
     * Example usage:
     *
          MSDataReader<CachedMzMLConsumer, Peak1D, ChromatogramPeak> exp_reader;
          MzMLFile mz_data_file;
          exp_reader.setConsumer(cacher);
          mz_data_file.load(in, exp_reader);
     *
    */
    template <typename ConsumerT, typename PeakT = Peak1D, typename ChromatogramPeakT = ChromatogramPeak>
    class OPENMS_DLLAPI MSDataReader :
      public MSDataIOInterfaceBase<PeakT, ChromatogramPeakT>
    {

  public:

      /// Constructor
      MSDataReader() {}

      // Need implementation, but we do not care about them
      inline void reserve(Size /* s */) {}
      inline void reserveSpaceChromatograms(Size /* s */) {}
      inline void reserveSpaceSpectra(Size /* s */) {}

      /// Resets all internal values
      void reset()
      {
        ExperimentalSettings::operator=(ExperimentalSettings()); //reset meta info
      }

      /// adds a spectrum to the consumer and keeps the meta-data (SpectrumSettings)
      void addSpectrum(const MSSpectrum<PeakT> & spectrum)
      {
        consumer->consumeSpectrum(spectrum);

        // We copy the meta-data of the spectrum
        MSSpectrum<PeakT> cpy = spectrum;
        cpy.clear(false);
        this->spectra_.push_back(cpy);
      }

      /// adds a chromatogram to the consumer and keeps the meta-data (ChromatogramSettings)
      void addChromatogram(const MSChromatogram<ChromatogramPeakT> & chromatogram)
      {
        consumer->consumeChromatogram(chromatogram);

        // We copy the meta-data of the chromatogram
        MSChromatogram<ChromatogramPeakT> cpy = chromatogram;
        cpy.clear(false);
        this->chromatograms_.push_back(cpy);
      }

      /// returns the list with the spectra settings
      const std::vector<MSSpectrum<PeakT> > & getSpectraSettings() const
      {
        return this->spectra_;
      }

      /// returns the list with the chromatogram settings
      const std::vector<MSChromatogram<ChromatogramPeakT> > & getChromatogramSettings() const
      {
        return this->chromatograms_;
      }

      inline void setConsumer(boost::shared_ptr<ConsumerT> c) { consumer = c; }

    protected:
      boost::shared_ptr<ConsumerT> consumer;
    };

    /*
     * Can count how many spectra and chromatograms are present in an mzML file. 
     *
     * This class relies on the fact that each spectrum and chromatogram will be
     * added to a writer through the addSpectrum and addChromatogram calls. The
     * drawback of this method (compared to the MSDataCounterReserve) is
     * that it is much slower. 
     *
     * Example usage:
     *
          MzMLFile f;
          MSDataCounter<Peak1D> exp_cnt;
          f.load(in, exp_cnt);
          // Result contained in exp_cnt.spectraCounts and exp_cnt.chromatogramCounts
     *
    */
    template <typename PeakT = Peak1D, typename ChromatogramPeakT = ChromatogramPeak>
    class OPENMS_DLLAPI MSDataCounter :
      public MSDataIOInterfaceBase<PeakT, ChromatogramPeakT>
    {

  public:

      /// Constructor
      MSDataCounter():
          spectraCounts(0),
          chromatogramCounts(0)
      {}

      inline void reserve(Size /* s */) {}
      inline void reserveSpaceChromatograms(Size /* s */) {}
      inline void reserveSpaceSpectra(Size /* s */) {}
      void reset() {}

      void addSpectrum(const MSSpectrum<PeakT> & /* spectrum */) { spectraCounts++; }
      void addChromatogram(const MSChromatogram<ChromatogramPeakT> & /* chromatogram */) { chromatogramCounts++; }

      Size spectraCounts;
      Size chromatogramCounts;
    };

    /*
     * Can count how many spectra and chromatograms are present in an mzML file. 
     *
     * This class relies on the fact that the spectrumList and chromatogramList
     * count attributes are accurate and that the MzMLHandler will try to reserve
     * appropriate space for them. 
     *
     * Example usage:
     *
          MzMLFile f;
          MSDataCounterReserve<Peak1D> exp_cnt;
          f.getOptions().addMSLevel(-1);
          f.load(in, exp_cnt);
          // Result contained in exp_cnt.spectraCounts and exp_cnt.chromatogramCounts
     *
    */
    template <typename PeakT = Peak1D, typename ChromatogramPeakT = ChromatogramPeak>
    class OPENMS_DLLAPI MSDataCounterReserve :
      public MSDataIOInterfaceBase<PeakT, ChromatogramPeakT>
    {

  public:

      /// Constructor
      MSDataCounterReserve():
          spectraCounts(0),
          chromatogramCounts(0)
      {}

      // grab the size of the chromatogram/spectra vector from the reserve calls
      inline void reserve(Size s) { spectraCounts = s; }
      inline void reserveSpaceSpectra(Size s) { spectraCounts = s; }
      inline void reserveSpaceChromatograms(Size s) { chromatogramCounts = s; }

      void reset() {}
      void addSpectrum(const MSSpectrum<PeakT> & /* spectrum */) {}
      void addChromatogram(const MSChromatogram<ChromatogramPeakT> & /* chromatogram */) {}

      Size spectraCounts;
      Size chromatogramCounts;
    };

    /*
     * Is able to consume calls to consumeSpectrum and consumeChromatogram.
     *
     * It will write spectra and chromatogram data to disk (location defined in
     * filename) in a binary format. Please create with a filename and then open
     * a file with the number of experiments/chromatograms.
     *
    */
    class OPENMS_DLLAPI CachedMzMLConsumer :
      public CachedmzML 
    {
      typedef MSSpectrum<Peak1D> SpectrumType;
      typedef MSChromatogram<ChromatogramPeak> ChromatogramType;

    public:
      /// Default constructor
      CachedMzMLConsumer(String filename) :
        ofs(filename.c_str(), std::ios::binary),
        spectra_written(0),
        chromatograms_written(0),
        spectra_expected(0),
        chromatograms_expected(0)
      {
      }

      /// Default destructor
      ~CachedMzMLConsumer()
      {
        ofs.close();
      }

      /// Write complete spectra as a dump to the disk
      void openFile(Size exp_size, Size chrom_size)
      {
        spectra_expected = exp_size;
        chromatograms_expected = chrom_size;

        int magic_number = MAGIC_NUMBER;
        ofs.write((char*)&magic_number, sizeof(magic_number));
        ofs.write((char*)&exp_size, sizeof(exp_size));
        ofs.write((char*)&chrom_size, sizeof(chrom_size));
      }

      void consumeSpectrum(const SpectrumType & s)
      {
        if (spectra_written >= spectra_expected || chromatograms_written > 0)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                  "Cannot write spectra, reached expected spectra or have already written chromagrams.");
        }
        writeSpectrum_(s, ofs);
        spectra_written++;
      }

      void consumeChromatogram(const ChromatogramType & c)
      {
        if (chromatograms_written >= chromatograms_expected || spectra_written != spectra_expected)
        {
          throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                  "Cannot write spectra, reached expected spectra or have already written chromagrams.");
        }
        writeChromatogram_(c, ofs);
        chromatograms_written++;
      }

    protected:
      std::ofstream ofs;
      Size spectra_written;
      Size chromatograms_written;
      Size spectra_expected;
      Size chromatograms_expected;

    };

    void noop (MSSpectrum<Peak1D> & /* s */) {;}

    template <typename MapType>
    class OPENMS_DLLAPI MzMLConsumer :
      public Internal::MzMLHandler<MapType>
    {

    public:
      typedef typename MapType::SpectrumType SpectrumType;
      typedef typename MapType::ChromatogramType ChromatogramType;

      /// Default constructor
      MzMLConsumer(String filename, ProgressLogger logger) :
        Internal::MzMLHandler<MapType>(MapType(), filename, MzMLFile().getVersion(), logger),
        ofs(filename.c_str()), 
        started_writing(false),
        spectra_written(0),
        chromatograms_written(0),
        spectra_expected(-1),
        chromatograms_expected(-1),
        add_dataprocessing_(false)
      {
        validator_ = new Internal::MzMLValidator(this->mapping_, this->cv_);
        sproptr_ = &noop; // setting default processing action to noop
      }

      /// Default destructor
      ~MzMLConsumer()
      {
        //--------------------------------------------------------------------------------------------
        //cleanup
        //--------------------------------------------------------------------------------------------
        ofs << "\t\t</spectrumList>\n";
        Internal::MzMLHandler<MapType>::writeFooter_(ofs);

        delete validator_;

        ofs.close();
      }

      /// Write complete spectra as a dump to the disk
      void setExperimentalSettings(ExperimentalSettings& exp)
      {
        settings = exp;
      }

      void setExpectedSize(Size expectedSpectra, Size expectedChromatograms)
      {
        spectra_expected = expectedSpectra;
        chromatograms_expected = expectedChromatograms;
      }

      void consumeSpectrum(const SpectrumType & s)
      {

        // Create copy and add dataprocessing if required
        SpectrumType scpy = s;
        processSpectrum_(scpy);
        if (add_dataprocessing_)
        {
          scpy.getDataProcessing().push_back(additional_dataprocessing_);
        }

        if (!started_writing)
        {
          // this is the first spectrum -> start writing the header
          // We also need to modify the map and add this dummy spectrum in
          // order to write the header correctly
          MapType dummy;
          dummy = settings;
          dummy.addSpectrum(scpy);

          //--------------------------------------------------------------------
          //header
          //--------------------------------------------------------------------
          writeHeader_(ofs, dummy, dps, *validator_);
          ofs << "\t\t<spectrumList count=\"" << spectra_expected << "\" defaultDataProcessingRef=\"dp_sp_0\">\n";
          started_writing = true;
        }
        bool renew_native_ids = false;
        Internal::MzMLHandler<MapType>::writeSpectrum_(ofs, scpy,
                spectra_written++, *validator_, renew_native_ids, dps);
      }

      void setSpectraProcessingPtr( void (*sproptr)(SpectrumType&) )
      {
        sproptr_ = sproptr;
      }

      void addDataProcessing(DataProcessing d)
      {
        additional_dataprocessing_ = d;
        add_dataprocessing_ = true;
      }

      void consumeChromatogram(const ChromatogramType & /* c */)
      {
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
    protected:

      virtual void processSpectrum_(SpectrumType & s)
      {
        // apply a function to it before writing ... 
        (*sproptr_)(s);
      }

      std::ofstream ofs;
      bool started_writing;
      Size spectra_written;
      Size chromatograms_written;
      Size spectra_expected;
      Size chromatograms_expected;
      bool add_dataprocessing_;

      ExperimentalSettings settings;
      Internal::MzMLValidator * validator_;
      std::vector<std::vector<DataProcessing> > dps;
      DataProcessing additional_dataprocessing_;

      void (*sproptr_)(SpectrumType&);

    };

  } // Internal

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
    OpenSwath::SpectrumAccessPtr cacheFile(String in, String tmp, String tmp_fname)
    {
      cacheFile_(in, tmp, tmp_fname);

      // Now re-load the cached file at tmp + tmp_fname and get a
      // SpectrumAccessPtr to return to caller.
      MSExperiment<Peak1D> exp;
      MzMLFile mzmlfile;
      mzmlfile.load(tmp + tmp_fname, exp);
      return SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(exp);
    }

  protected:

    /**
        @brief Reads and input mzML file and caches it to a temporary location.
    */
    void cacheFile_(String in, String tmp, String tmp_fname)
    {
      String cached_file = tmp + tmp_fname + ".cached";
      String meta_file = tmp + tmp_fname;

      // Load the file only to count the number of spectra/chromatograms
      InternalCaching::MSDataCounterReserve<Peak1D> exp_cnt;
      {
        MzMLFile f;
        f.getOptions().addMSLevel(-1);
        f.load(in, exp_cnt);
      }

      // Open a cacher file with the known amounts of spectra/chromatograms at
      // the temporary location "cached_file"
      boost::shared_ptr<InternalCaching::CachedMzMLConsumer> cacher(new InternalCaching::CachedMzMLConsumer(cached_file));
      cacher->openFile(exp_cnt.spectraCounts, exp_cnt.chromatogramCounts);

      // Create a new file acceptor
      InternalCaching::MSDataReader<InternalCaching::CachedMzMLConsumer, Peak1D, ChromatogramPeak> exp_reader;
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
