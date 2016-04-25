#include "ACTS/Plugins/DD4hepConverters/DetExtension.h"

Acts::DetExtension::DetExtension() :
Acts::IDetExtension(),
m_segmentation(nullptr),
m_shape(Acts::ShapeType::None),
m_supportStructure(),
m_modules()
{}

Acts::DetExtension::~DetExtension()
{}