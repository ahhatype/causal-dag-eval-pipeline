"""Minimal reader for the GRKM container (format version 1: raw sections)."""
import struct, zlib

def sections(buf):
    assert buf[:4] == b"GRKM", "not a model file"
    ver, endian, count = struct.unpack_from("<IB3xI", buf, 4)
    assert endian == 1
    off = 16
    out = []
    for _ in range(count):
        tag = buf[off:off+16].rstrip(b"\0").decode()
        if ver == 1:
            (ln, crc) = struct.unpack_from("<QI", buf, off+16)
            start = off + 28
            payload = buf[start:start+ln]
            off = start + ln
        else:
            ln, crc, codec, raw = struct.unpack_from("<QIBQ", buf, off+16)
            start = off + 37
            payload = buf[start:start+ln]
            if codec != 0:
                raise NotImplementedError("compressed section %s (format 2)" % tag)
            off = start + ln
        assert zlib.crc32(payload) & 0xffffffff == crc, "CRC mismatch in " + tag
        out.append((tag, payload))
    return ver, out
