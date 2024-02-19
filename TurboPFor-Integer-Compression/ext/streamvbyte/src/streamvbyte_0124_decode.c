#include "streamvbyte.h"
#include "streamvbyte_isadetection.h"

#ifdef STREAMVBYTE_X64
#include "streamvbyte_shuffle_tables_0124_decode.h"
#endif

#include <string.h> // for memcpy



#ifdef STREAMVBYTE_X64 
STREAMVBYTE_TARGET_SSSE3
static inline __m128i _decode_avx(uint32_t key,
                                  const uint8_t *__restrict__ *dataPtrPtr) {
  uint8_t len;
  __m128i Data = _mm_loadu_si128((__m128i *)*dataPtrPtr);
  uint8_t *pshuf = (uint8_t *) &shuffleTable[key];
  __m128i Shuf = *(__m128i *)pshuf;
  len = lengthTable[key];
  Data = _mm_shuffle_epi8(Data, Shuf);
  *dataPtrPtr += len;
  return Data;
}
STREAMVBYTE_UNTARGET_REGION

STREAMVBYTE_TARGET_SSSE3
static inline void _write_avx(uint32_t *out, __m128i Vec) {
  _mm_storeu_si128((__m128i *)out, Vec);
}
STREAMVBYTE_UNTARGET_REGION

#endif // STREAMVBYTE_X64

static inline uint32_t _decode_data(const uint8_t **dataPtrPtr, uint8_t code) {
  const uint8_t *dataPtr = *dataPtrPtr;
  uint32_t val;

  if (code == 0) { // 0 byte
    val = 0;
  } else if (code == 1) { // 1 bytes
    val = (uint32_t)*dataPtr;
    dataPtr += 1;
  } else if (code == 2) { // 2 bytes
    val = 0;
    memcpy(&val, dataPtr, 2); // assumes little endian
    dataPtr += 2;
  } else { // code == 3, 4 bytes
    memcpy(&val, dataPtr, 4);
    dataPtr += 4;
  }

  *dataPtrPtr = dataPtr;
  return val;
}

static const uint8_t *svb_decode_scalar(uint32_t *outPtr, const uint8_t *keyPtr,
                                        const uint8_t *dataPtr,
                                        uint32_t count) {
  if (count == 0)
    return dataPtr; // no reads or writes if no data

  uint8_t shift = 0;
  uint32_t key = *keyPtr++;
  for (uint32_t c = 0; c < count; c++) {
    if (shift == 8) {
      shift = 0;
      key = *keyPtr++;
    }
    uint32_t val = _decode_data(&dataPtr, (key >> shift) & 0x3);
    *outPtr++ = val;
    shift += 2;
  }

  return dataPtr; // pointer to first unused byte after end
}

#ifdef STREAMVBYTE_X64
STREAMVBYTE_TARGET_SSSE3
static const uint8_t *svb_decode_avx_simple(uint32_t *out,
                                            const uint8_t *__restrict__ keyPtr,
                                     const uint8_t *__restrict__ dataPtr,
                                     uint64_t count) {

  uint64_t keybytes = count / 4; // number of key bytes
  __m128i Data;
  if (keybytes >= 8) {

    int64_t Offset = -(int64_t)keybytes / 8 + 1;

    const uint64_t *keyPtr64 = (const uint64_t *)keyPtr - Offset;
    uint64_t nextkeys;
    memcpy(&nextkeys, keyPtr64 + Offset, sizeof(nextkeys));
    for (; Offset != 0; ++Offset) {
      uint64_t keys = nextkeys;
      memcpy(&nextkeys, keyPtr64 + Offset + 1, sizeof(nextkeys));

      Data = _decode_avx((keys & 0xFF), &dataPtr);
      _write_avx(out, Data);
      Data = _decode_avx((keys & 0xFF00) >> 8, &dataPtr);
      _write_avx(out + 4, Data);

      keys >>= 16;
      Data = _decode_avx((keys & 0xFF), &dataPtr);
      _write_avx(out + 8, Data);
      Data = _decode_avx((keys & 0xFF00) >> 8, &dataPtr);
      _write_avx(out + 12, Data);

      keys >>= 16;
      Data = _decode_avx((keys & 0xFF), &dataPtr);
      _write_avx(out + 16, Data);
      Data = _decode_avx((keys & 0xFF00) >> 8, &dataPtr);
      _write_avx(out + 20, Data);

      keys >>= 16;
      Data = _decode_avx((keys & 0xFF), &dataPtr);
      _write_avx(out + 24, Data);
      Data = _decode_avx((keys & 0xFF00) >> 8, &dataPtr);
      _write_avx(out + 28, Data);

      out += 32;
    }
    {
      uint64_t keys = nextkeys;

      Data = _decode_avx((keys & 0xFF), &dataPtr);
      _write_avx(out, Data);
      Data = _decode_avx((keys & 0xFF00) >> 8, &dataPtr);
      _write_avx(out + 4, Data);

      keys >>= 16;
      Data = _decode_avx((keys & 0xFF), &dataPtr);
      _write_avx(out + 8, Data);
      Data = _decode_avx((keys & 0xFF00) >> 8, &dataPtr);
      _write_avx(out + 12, Data);

      keys >>= 16;
      Data = _decode_avx((keys & 0xFF), &dataPtr);
      _write_avx(out + 16, Data);
      Data = _decode_avx((keys & 0xFF00) >> 8, &dataPtr);
      _write_avx(out + 20, Data);

      keys >>= 16;
      Data = _decode_avx((keys & 0xFF), &dataPtr);
      _write_avx(out + 24, Data);
      Data = _decode_avx((keys & 0xFF00) >> 8, &dataPtr);
      _write_avx(out + 28, Data);

      out += 32;
    }
  }

  return dataPtr;
}
STREAMVBYTE_UNTARGET_REGION


#endif

// Read count 32-bit integers in maskedvbyte format from in, storing the result
// in out.  Returns the number of bytes read.
size_t streamvbyte_decode_0124(const uint8_t *in, uint32_t *out, uint32_t count) {
  if (count == 0)
    return 0;


  const uint8_t *keyPtr = in;               // full list of keys is next
  uint32_t keyLen = ((count + 3) / 4);      // 2-bits per key (rounded up)
  const uint8_t *dataPtr = keyPtr + keyLen; // data starts at end of keys

#ifdef STREAMVBYTE_X64
  if(streamvbyte_ssse3()) {
    dataPtr = svb_decode_avx_simple(out, keyPtr, dataPtr, count);
    out += count & ~ 31;
    keyPtr += (count/4) & ~ 7;
    count &= 31;
  }
#endif

  return svb_decode_scalar(out, keyPtr, dataPtr, count) - in;

}
