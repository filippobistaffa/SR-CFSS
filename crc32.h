#ifndef CRC32_H_
#define CRC32_H_

#include <stdlib.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C"
#endif
uint32_t crc32(const void *buf, size_t size);

#endif /* CRC32_H_ */
