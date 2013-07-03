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

#ifndef OPENMS_LOWMEM_MZMLCONSUMER
#define OPENMS_LOWMEM_MZMLCONSUMER

#include "OpenMS/ANALYSIS/OPENSWATH/CachedmzML.h"
#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>

/*
 * MzMLConsumer can consume calls to consumeSpectrum and consumeChromatogram.
 *
 * The idea is that MSDataReader will accept a consumer and forward all read
 * spectra/chromatograms to the consumer.
 *
*/

namespace OpenMS
{

  namespace Internal
  {
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

} // NS OpenMS

#endif
