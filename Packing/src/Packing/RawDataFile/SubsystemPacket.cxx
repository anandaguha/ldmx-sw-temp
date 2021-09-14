
#include "Packing/RawDataFile/SubsystemPacket.h"

#include "Packing/Utility/Mask.h"

namespace packing {
namespace rawdatafile {

void SubsystemPacket::read(Reader& r) {
  uint32_t word;
  r.read(&word, 1);

  subsys_id_ = (word >> 16) & mask<16>;
  uint32_t subsys_len = (word >> 1) & mask<15>;
  crc_ok_ = word & mask<1>;

  data_.reserve(subsys_len);
  r.read(data_.data(), subsys_len);

  r.read(&crc_, 1);
}

/*
void SubsystemPacket::write(std::ofstream& os) {
  static WordType word;

  word = ((subsys_id_ & mask<16>) << 16
        +(data_.size() & mask<15>) << 1
        +(crc_ok_ ? 1 : 0));

  word >> os;

  for (const auto& w : data_) w >> os;

  crc_ >> os;
}
*/

}  // namespace rawdatafile
}  // namespace packing
