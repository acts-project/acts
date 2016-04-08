#include "DD4hepDetectorElement/DetExtension.h"

Acts::DetExtension::DetExtension() :
Acts::IDetExtension(),
m_segmentation(nullptr),
m_shape(Acts::ShapeType::None)
{}

Acts::DetExtension::~DetExtension()
{}