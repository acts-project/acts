// Copyright (C) 2002-2019 CERN for the benefit of the ATLAS collaboration

/**
 * @file HTTHoughTransformTool.cxx
 * @author Riley Xu - riley.xu@cern.ch
 * @date October 31st, 2020
 * @brief See header file.
 */

#include "TrigHTTUtils/HTTTypes.h"
#include "TrigHTTConfig/IHTTEventSelectionSvc.h"
#include "TrigHTTObjects/HTTHit.h"
#include "TrigHTTMaps/ITrigHTTMappingSvc.h"
#include "TrigHTTMaps/HTTPlaneMap.h"
#include "TrigHTTMaps/HTTRegionMap.h"
#include "TrigHTTBanks/ITrigHTTBankSvc.h"
#include "TrigHTTBanks/HTTSectorBank.h"
#include "TrigHTTAlgorithms/HTTHoughTransformTool.h"

#include <sstream>
#include <cmath>
#include <algorithm>

static inline int quant(double min, double max, unsigned nSteps, double val);
static inline double unquant(double min, double max, unsigned nSteps, int step);
template <typename T>
static inline std::string to_string(std::vector<T> v);

//const double HTTHoughTransformTool::A = 0.0003;

///////////////////////////////////////////////////////////////////////////////
// AthAlgTool

HTTHoughTransformTool::HTTHoughTransformTool(const std::string& algname, const std::string &name, const IInterface *ifc) :
    AthAlgTool(algname, name, ifc),
    m_EvtSel("HTTEventSelectionSvc", name),
    m_HTTBankSvc("TrigHTTBankSvc", name),
    m_HTTMapping("TrigHTTMappingSvc", name)
{
    declareInterface<HTTRoadFinderToolI>(this);
    declareProperty("subRegion", m_subRegion, "Sub region of this transform, or -1 for full region");
    declareProperty("phi_min", m_parMin.phi, "Minimum phi of transform");
    declareProperty("phi_max", m_parMax.phi, "Maximum phi of transform");
    declareProperty("qpT_min", m_parMin.qOverPt, "Min q/pT of transform");
    declareProperty("qpT_max", m_parMax.qOverPt, "Max q/pT of transform");
    declareProperty("d0_min", m_parMin.d0, "Min d0 of transform");
    declareProperty("d0_max", m_parMax.d0, "Max d0 of transform");
    declareProperty("nBins_x", m_imageSize_x);
    declareProperty("nBins_y", m_imageSize_y);
    declareProperty("convolution", m_conv, "Vector containing the convolution filter");
    declareProperty("combine_layers", m_combineLayers, "Vector containing the layers we want to combine");
    declareProperty("scale", m_binScale, "Vector containing the scales for each layers");
    declareProperty("convSize_x", m_convSize_x);
    declareProperty("convSize_y", m_convSize_y);
    declareProperty("hitExtend_x", m_hitExtend_x, "Hit lines will fill extra bins in x by this amount on each side, size == nLayers");
    declareProperty("threshold", m_threshold, "Minimum value post-convolution to accept as a road (inclusive)");
    declareProperty("traceHits", m_traceHits, "Trace each hit that goes in a bin. Disabling this will save memory/time since each bin doesn't have to store all its hits but the roads created won't have hits from convolution, etc.");
    declareProperty("localMaxWindowSize", m_localMaxWindowSize, "Only create roads that are a local maximum within this window size. Set this to 0 to turn off local max filtering");
    declareProperty("useSectors", m_useSectors, "Will reverse calculate the sector for track-fitting purposes");
    declareProperty("IdealGeoRoads", m_idealGeoRoads, "Set sectors to use ideal geometry fit constants");
    declareProperty("fieldCorrection", m_fieldCorrection, "Apply corrections to hough equation due to field nonuniformity");
}


StatusCode HTTHoughTransformTool::initialize()
{
    // Debug
    ATH_MSG_INFO("Image size: " << m_imageSize_x << " x " << m_imageSize_y);
    ATH_MSG_INFO("Convolution size: " << m_convSize_x << " x " << m_convSize_y);
    ATH_MSG_INFO("Convolution: " << to_string(m_conv));
    ATH_MSG_INFO("Hit Extend: " << to_string(m_hitExtend_x));

    // Retrieve info
    ATH_CHECK(m_HTTBankSvc.retrieve());
    ATH_CHECK(m_HTTMapping.retrieve());
    m_nLayers = m_HTTMapping->PlaneMap_1st()->getNLogiLayers();

    // Error checking
    // TODO check bounds are set correctly
    bool ok = false;
    if (!m_imageSize_x || !m_imageSize_y)
        ATH_MSG_FATAL("initialize() Image size must be greater than 0");
    else if (m_conv.size() != m_convSize_x * m_convSize_y)
        ATH_MSG_FATAL("initialize() Convolution sizes don't match");
    else if (!m_conv.empty() && (m_convSize_x % 2 == 0 || m_convSize_y % 2 == 0))
        ATH_MSG_FATAL("initialize() Convolution sizes must be odd");
    else if (m_hitExtend_x.size() % m_nLayers)
        ATH_MSG_FATAL("initialize() Hit extentsion list must have size % nLayers");
    else if (!m_combineLayers.empty() && m_combineLayers.size() != m_nLayers)
        ATH_MSG_FATAL("initialize() Combine layers list must have size = nLayers");
    else if (m_threshold.size() % 2 != 1)
        ATH_MSG_FATAL("initialize() Threshold size must be odd");
    else if (!m_binScale.empty() && m_binScale.size() != m_nLayers)
        ATH_MSG_FATAL("initialize() Bin scale list must have size = nLayers");
    else if (std::any_of(m_binScale.begin(), m_binScale.end(), [&](unsigned i){ return m_imageSize_y % i != 0; }))
        ATH_MSG_FATAL("initialize() The imagesize is not divisible by scale");
    else
        ok = true;
    if (!ok) return StatusCode::FAILURE;

    // Warnings / corrections
    if (m_localMaxWindowSize && !m_traceHits)
    {
        ATH_MSG_WARNING("initialize() localMaxWindowSize requires tracing hits, turning on automatically");
        m_traceHits = true;
    }
    if (m_idealGeoRoads)
    {
        if (m_useSectors)
        {
            ATH_MSG_WARNING("initialize() idealGeoRoads conflicts with useSectors, switching off HTT sector matching");
            m_useSectors = false;
        }
        if (!m_traceHits)
        {
            ATH_MSG_WARNING("initialize() idealGeoRoads requires tracing hits, turning on automatically");
            m_traceHits = true;
        }
    }
    if (m_binScale.empty()) m_binScale.resize(m_nLayers, 1);

    // Fill convenience variables
    m_step_x = (m_parMax[m_par_x] - m_parMin[m_par_x]) / m_imageSize_x;
    m_step_y = (m_parMax[m_par_y] - m_parMin[m_par_y]) / m_imageSize_y;
    for (unsigned i = 0; i <= m_imageSize_x; i++)
        m_bins_x.push_back(unquant(m_parMin[m_par_x], m_parMax[m_par_x], m_imageSize_x, i));
    for (unsigned i = 0; i <= m_imageSize_y; i++)
        m_bins_y.push_back(unquant(m_parMin[m_par_y], m_parMax[m_par_y], m_imageSize_y, i));

    // Initialize combine layers
    if (!m_combineLayers.empty())
    {
        m_nCombineLayers = *std::max_element(m_combineLayers.begin(), m_combineLayers.end()) + 1;
        m_combineLayer2D.resize(m_nCombineLayers);
        for (unsigned i = 0; i < m_combineLayers.size(); i++)
            m_combineLayer2D[m_combineLayers[i]].push_back(i);
    }
    else
    {
        m_nCombineLayers = m_nLayers;
        for (unsigned i = 0; i < m_nLayers; i++)
            m_combineLayer2D.push_back({ i });
    }

    return StatusCode::SUCCESS;
}


StatusCode HTTHoughTransformTool::finalize()
{
    return StatusCode::SUCCESS;
}


///////////////////////////////////////////////////////////////////////////////
// Main Algorithm

StatusCode HTTHoughTransformTool::getRoads(const std::vector<const HTTHit*> & hits, std::vector<HTTRoad*> & roads) 
{
    roads.clear();
    m_roads.clear();

    m_image = createImage(hits);
    if (!m_conv.empty()) m_image = convolute(m_image);

    for (unsigned y = 0; y < m_imageSize_y; y++)
        for (unsigned x = 0; x < m_imageSize_x; x++)
            if (passThreshold(m_image, x, y))
            {
                if (m_traceHits)
                    addRoad(m_image(y, x).second, x, y);
                else
                    addRoad(hits, x, y);
            }

    roads.reserve(m_roads.size());
    for (HTTRoad_Hough & r : m_roads) roads.push_back(&r);

    return StatusCode::SUCCESS;
}

HTTHoughTransformTool::Image HTTHoughTransformTool::createLayerImage(std::vector<unsigned> const & layers, std::vector<HTTHit const *> const & hits, unsigned const scale) const
{
    Image image(m_imageSize_y, m_imageSize_x);

    for (HTTHit const * hit : hits)
    {
        if (std::find(layers.begin(), layers.end(), hit->getLayer()) == layers.end()) continue;
        if (m_subRegion >= 0 && !m_HTTMapping->SubRegionMap()->isInRegion(m_subRegion, *hit)) continue;

        // This scans over y (pT) because that is more efficient in memory, in C.
        // Unknown if firmware will want to scan over x instead.
        unsigned new_size_y  = m_imageSize_y / scale;
        for (unsigned y_ = 0; y_ < new_size_y; y_++)
        {
            unsigned y_bin_min = scale * y_;
            unsigned y_bin_max = scale * (y_ + 1);

            // Find the min/max x bins
            auto xBins = yToXBins(y_bin_min, y_bin_max, hit);

            // Update the image
            for (unsigned y = y_bin_min; y < y_bin_max; y++)
                for (unsigned x = xBins.first; x < xBins.second; x++)
                {
                    image(y, x).first++;
                    if (m_traceHits) image(y, x).second.insert(hit);
                }
        }
    }

    return image;
}

HTTHoughTransformTool::Image HTTHoughTransformTool::createImage(std::vector<HTTHit const *> const & hits) const
{
    Image image(m_imageSize_y, m_imageSize_x);

    for (unsigned i = 0; i < m_nCombineLayers; i++)
    {
        Image layerImage = createLayerImage(m_combineLayer2D[i], hits, m_binScale[i]);
        for (unsigned x = 0; x < m_imageSize_x; ++x)
            for (unsigned y = 0; y < m_imageSize_y; ++y)
                if (layerImage(y, x).first > 0)
                {
                    image(y, x).first++;
                    image(y, x).second.insert(layerImage(y, x).second.begin(), layerImage(y, x).second.end());
                }
    }

    return image;
}

HTTHoughTransformTool::Image HTTHoughTransformTool::convolute(Image const & image) const
{
    Image out(m_imageSize_y, m_imageSize_x);

    for (unsigned y0 = 0; y0 < m_imageSize_y; y0++)     // Loop over out
        for (unsigned x0 = 0; x0 < m_imageSize_x; x0++) //
            for (unsigned r = 0; r < m_convSize_y; r++)     // Loop over conv
                for (unsigned c = 0; c < m_convSize_x; c++) //
                {
                    int y = -static_cast<int>(m_convSize_y) / 2 + r + y0; // Indices of input
                    int x = -static_cast<int>(m_convSize_x) / 2 + c + x0; //

                    if (y >= 0 && y < static_cast<int>(m_imageSize_y) && x >= 0 && x < static_cast<int>(m_imageSize_x))
                    {
                        int val = m_conv[r * m_convSize_x + c] * image(y, x).first;
                        if (val > 0)
                        {
                            out(y0, x0).first += val;
                            out(y0, x0).second.insert(image(y, x).second.begin(), image(y, x).second.end());
                        }
                    }
                }

    return out;
}

bool HTTHoughTransformTool::passThreshold(Image const & image, unsigned x, unsigned y) const
{
    // Pass window threshold
    unsigned width = m_threshold.size() / 2;
    if (x < width || (image.size(1) - x) < width) return false;
    for (unsigned i = 0; i < m_threshold.size(); i++)
        if (image(y, x - width + i).first < m_threshold[i]) return false;

    // Pass local-maximum check
    if (m_localMaxWindowSize)
        for (int j = -m_localMaxWindowSize; j <= m_localMaxWindowSize; j++)
            for (int i = -m_localMaxWindowSize; i <= m_localMaxWindowSize; i++)
            {
                if (i == 0 && j == 0) continue;
                if (y + j < image.size(0) && x + i < image.size(1))
                {
                    if (image(y+j, x+i).first > image(y, x).first) return false;
                    if (image(y+j, x+i).first == image(y, x).first)
                    {
                        if (image(y+j, x+i).second.size() > image(y, x).second.size()) return false;
                        if (image(y+j, x+i).second.size() == image(y, x).second.size()
                            && j <= 0 && i <= 0) return false; // favor bottom-left (low phi, low neg q/pt)
                    }
                }
            }

    return true;
}

///////////////////////////////////////////////////////////////////////////////
// Helpers


// Quantizes val, given a range [min, max) split into nSteps. Returns the bin below.
static inline int quant(double min, double max, unsigned nSteps, double val)
{
    return static_cast<int>((val - min) / (max - min) * nSteps);
}

// Returns the lower bound of the bin specified by step
static inline double unquant(double min, double max, unsigned nSteps, int step)
{
    return min + (max - min) * step / nSteps;
}

template <typename T>
static inline std::string to_string(std::vector<T> v)
{
    std::ostringstream oss;
    oss << "[";
    if (!v.empty())
    {
        std::copy(v.begin(), v.end()-1, std::ostream_iterator<T>(oss, ", "));
        oss << v.back();
    }
    oss << "]";
    return oss.str();
}


double HTTHoughTransformTool::fieldCorrection(unsigned region, double qpt, double r)
{
    r = r / 1000; // convert to meters
    if (region == 3)
    {
        double cor = 0.1216 * r * r - 0.0533 * r + 0.0069;
        return -cor * qpt;
    }
    else if (region == 4)
    {
        double cor = 0.4265 * r * r - 0.0662 * r + 0.0036;
        return -cor * qpt;
    }
    else return 0;
}

double HTTHoughTransformTool::yToX(double y, HTTHit const * hit) const
{
    double x = 0;

    if (m_par_x == HTTTrackPars::IPHI && m_par_y == HTTTrackPars::IHIP)
    {
        double r = hit->getR(); // mm
        double phi_hit = hit->getGPhi(); // radians
        double d0 = std::isnan(m_parMin.d0) ? 0 : m_parMin.d0; // mm, assume min = max
        x = asin(r * A * y - d0 / r) + phi_hit;

        if (m_fieldCorrection) x += fieldCorrection(m_EvtSel->getRegionID(), y, r);
    }
    else
    {
        ATH_MSG_ERROR("yToX() not defined for the current m_par selection");
    }

    return x;
}

// Find the min/max x bins of the hit's line, in each y bin. Max is exclusive.
// Note this assumes yToX is monotonic. Returns {0, 0} if hit lies out of bounds.
std::pair<unsigned, unsigned> HTTHoughTransformTool::yToXBins(size_t yBin_min, size_t yBin_max, HTTHit const * hit) const
{
    // Get float values
    double x_min = yToX(m_bins_y[yBin_min], hit);
    double x_max = yToX(m_bins_y[yBin_max], hit);
    if (x_min > x_max) std::swap(x_min, x_max);
    if (x_max < m_parMin[m_par_x] || x_min > m_parMax[m_par_x])
        return { 0, 0 }; // out of bounds

    // Get bins
    int x_bin_min = quant(m_parMin[m_par_x], m_parMax[m_par_x], m_imageSize_x, x_min);
    int x_bin_max = quant(m_parMin[m_par_x], m_parMax[m_par_x], m_imageSize_x, x_max) + 1; // exclusive

    // Extend bins
    unsigned extend = getExtension(yBin_min, hit->getLayer());
    x_bin_min -= extend;
    x_bin_max += extend;

    // Clamp bins
    if (x_bin_min < 0) x_bin_min = 0;
    if (x_bin_max > static_cast<int>(m_imageSize_x)) x_bin_max = m_imageSize_x;

    return { x_bin_min, x_bin_max };
}

// We allow variable extension based on the size of m_hitExtend_x. See comments below.
unsigned HTTHoughTransformTool::getExtension(unsigned y, unsigned layer) const
{
    if (m_hitExtend_x.size() == m_nLayers) return m_hitExtend_x[layer];
    if (m_hitExtend_x.size() == m_nLayers * 2)
    {
        // different extension for low pt vs high pt, split in half but irrespective of sign
        // first nLayers entries of m_hitExtend_x is for low pt half, rest are for high pt half
        if (y < m_imageSize_y / 4 || y > 3 * m_imageSize_y / 4) return m_hitExtend_x[layer];
        return m_hitExtend_x[m_nLayers + layer];
    }
    return 0;
}

void HTTHoughTransformTool::matchIdealGeoSector(HTTRoad_Hough & r) const
{
    float pt = r.getY()*0.001; // convert to MeV
    auto bounds = std::equal_range(htt::QOVERPT_BINS.begin(),htt::QOVERPT_BINS.end(),pt);
    int sectorbin = bounds.first-htt::QOVERPT_BINS.begin()-1;

    // those bins are for tracks between the values, can't be below first value or more than the last value
    if (sectorbin < 0) sectorbin = 0;
    if (sectorbin > static_cast<int>(htt::QOVERPT_BINS.size()-2)) sectorbin =  htt::QOVERPT_BINS.size()-2;
    std::vector<module_t> modules;

    for (unsigned int il = 0; il < r.getNLayers(); il++) {
        if (r.getNHits_layer()[il] == 0) {
            modules.push_back(-1);
            layer_bitmask_t wc_layers = r.getWCLayers();
            wc_layers |= (0x1 << il);
            r.setWCLayers(wc_layers);

            HTTHit *wcHit = new HTTHit();
            wcHit->setHitType(HitType::wildcard);
            wcHit->setLayer(il);
            wcHit->setDetType(m_HTTMapping->PlaneMap_1st()->getDetType(il));
            std::vector<const HTTHit*> wcHits;
            wcHits.push_back(wcHit);
            r.setHits(il,wcHits);
        }
        else
            modules.push_back(sectorbin);
    }
    HTTSectorBank* sectorbank = m_HTTBankSvc->SectorBank_1st();
    r.setSector(sectorbank->findSector(modules));
}

// Creates a road from hits that pass through the given bin (x, y), and pushes it onto m_roads
void HTTHoughTransformTool::addRoad(std::vector<std::vector<const HTTHit*>> const & hits, layer_bitmask_t hitLayers, unsigned x, unsigned y)
{
    m_roads.emplace_back();
    HTTRoad_Hough & r = m_roads.back();

    r.setRoadID(m_roads.size() - 1);
    r.setPID(y * m_imageSize_y + x);
    r.setHits(hits);
    if (m_useSectors) r.setSector(m_HTTBankSvc->SectorBank_1st()->findSector(hits));
    else if (m_idealGeoRoads) matchIdealGeoSector(r);
    r.setHitLayers(hitLayers);
    r.setSubRegion(m_subRegion);
    r.setX(m_bins_x[x] + m_step_x/2);
    r.setY(m_bins_y[y] + m_step_y/2);
    r.setXBin(x);
    r.setYBin(y);
}


// Creates a road from hits that pass through the given bin (x, y), and pushes it onto m_roads
void HTTHoughTransformTool::addRoad(std::unordered_set<const HTTHit*> const & hits, unsigned x, unsigned y)
{
    layer_bitmask_t hitLayers = 0;
    for (HTTHit const * hit : hits)
        hitLayers |= 1 << hit->getLayer();

    auto sorted_hits = ::sortByLayer(hits);
    sorted_hits.resize(m_nLayers); // If no hits in last layer, return from sortByLayer will be too short

    addRoad(sorted_hits, hitLayers, x, y);
}

// Use this version of addRoad when hit tracing is turned off
void HTTHoughTransformTool::addRoad(std::vector<const HTTHit*> const & hits, unsigned x, unsigned y)
{
    // Get the road hits
    std::vector<HTTHit const *> road_hits;
    layer_bitmask_t hitLayers = 0;
    for (const HTTHit * hit : hits)
    {
        if (m_subRegion >= 0 && !m_HTTMapping->SubRegionMap()->isInRegion(m_subRegion, *hit)) continue;

        // Find the min/max y bins (after scaling)
        unsigned y_bin_min = (y / m_binScale[hit->getLayer()]) * m_binScale[hit->getLayer()];
        unsigned y_bin_max = y_bin_min + m_binScale[hit->getLayer()];

        // Find the min/max x bins
        auto xBins = yToXBins(y_bin_min, y_bin_max, hit);
        if (x >= xBins.first && x < xBins.second)
        {
            road_hits.push_back(hit);
            hitLayers |= 1 << hit->getLayer();
        }
    }

    auto sorted_hits = ::sortByLayer(road_hits);
    sorted_hits.resize(m_nLayers); // If no hits in last layer, return from sortByLayer will be too short

    addRoad(sorted_hits, hitLayers, x, y);
}

