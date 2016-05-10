#include "ACTS/Plugins/DD4hepPlugins/DetExtension.h"

Acts::DetExtension::DetExtension() :
Acts::IDetExtension(),
m_segmentation(nullptr),
m_shape(Acts::ShapeType::None),
m_supportStructure(),
m_modules()
{}

Acts::DetExtension::DetExtension(const DD4hep::Geometry::Segmentation segmentation) :
Acts::IDetExtension(),
m_segmentation(segmentation),
m_shape(Acts::ShapeType::None),
m_supportStructure(),
m_modules()
{}

Acts::DetExtension::DetExtension(ShapeType shape) :
Acts::IDetExtension(),
m_segmentation(nullptr),
m_shape(shape),
m_supportStructure(),
m_modules()
{}

Acts::DetExtension::DetExtension(const DD4hep::Geometry::DetElement support) :
Acts::IDetExtension(),
m_segmentation(nullptr),
m_shape(Acts::ShapeType::None),
m_supportStructure(support),
m_modules()
{}

Acts::DetExtension::DetExtension(std::vector<Module> mod) :
Acts::IDetExtension(),
m_segmentation(nullptr),
m_shape(Acts::ShapeType::None),
m_supportStructure(),
m_modules(mod)
{}
Acts::DetExtension::DetExtension(const DD4hep::Geometry::DetElement support, std::vector<Module> mod) :
Acts::IDetExtension(),
m_segmentation(nullptr),
m_shape(Acts::ShapeType::None),
m_supportStructure(support),
m_modules(mod)
{}


Acts::DetExtension::~DetExtension()
{}