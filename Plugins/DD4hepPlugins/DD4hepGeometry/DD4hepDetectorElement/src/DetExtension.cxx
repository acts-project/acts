#include "DD4hepDetectorElement/DetExtension.h"

Add4hep::DetExtension::DetExtension() :
Add4hep::IDetExtension(),
m_shape(Add4hep::ShapeType::None),
m_segmentation(nullptr)
{}

Add4hep::DetExtension::~DetExtension()
{}