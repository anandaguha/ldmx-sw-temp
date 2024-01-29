#pragma once

#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Definitions/Algebra.hpp"


namespace tracking::geo {

/**
 * The geometry identifier will return vol=0 and lay=0 when it is not valid.
 * The geometry identifier is only valid after building the tracking geometry.
 */

inline unsigned int unpackGeometryIdentifier(const Acts::GeometryIdentifier& geoId) {
  
  unsigned int volumeId = geoId.volume();
  unsigned int layerId  = geoId.layer() / 2;
  unsigned int sensorId = geoId.sensitive() - 1;
  unsigned int surfaceId = volumeId * 1000 + layerId * 100 + sensorId;
  
  return surfaceId;
}

// I use the same ATLAS convention (opposite to MPII)
//deltaR = (ru, rv, rw)
//          /  1    -rw   rv  \
//deltaR => |  rw    1   -ru  |
//          \ -rv    ru   1   /

inline Acts::RotationMatrix3 deltaRot(Acts::Vector3 deltaR) {

  //This is fine because RotationMatrix3 doesn't need to be symmetric
  Acts::RotationMatrix3 rot = Acts::RotationMatrix3::Identity();
  rot(0,1) = -deltaR(2);
  rot(0,2) =  deltaR(1);
  rot(1,2) = -deltaR(0);
  
  rot(1,0) = -rot(0,1);
  rot(2,0) = -rot(0,2);
  rot(2,1) = -rot(1,2);
  
  return rot;
  
}

}
//Acts::GeometryIdentifier packGeometryIdentifier(const unsigned int surfaceId) {}







