#ifndef __CHARB_HPP__
#define __CHARB_HPP__

#include <cstdlib>
#include <cstdio>
#include <cstdarg>
#include <stdexcept>
#include <string>
#include <cstring>

class charb {
  char *_base;
  char *_end;
  char *_ptr;

public:
  charb() : _base(0), _end(0), _ptr(0) { }
  explicit charb(size_t s) : _base(0), _end(0), _ptr(0) {
    ensure(s);
    *_base = '\0';
  }
  charb(const charb &rhs) : _base(0), _end(0), _ptr(0) {
    ensure(rhs.capacity());
    memcpy(_base, rhs._base, capacity());
  }
  charb(const char *str) : _base(0), _end(0), _ptr(0) {
    size_t len = strlen(str);
    ensure(len + 1);
    memcpy(_base, str, len + 1);
    _ptr = _base + len;
  }
  charb(const char *str, size_t len) : _base(0), _end(0), _ptr(0) {
    ensure(len + 1);
    memcpy(_base, str, len + 1);
    _ptr = _base + len;
  }
  charb(const std::string &s) : _base(0), _end(0), _ptr(0) {
    ensure(s.size() + 1);
    memcpy(_base, s.c_str(), s.size() + 1);
    _ptr = _base + s.size();
  }
  virtual ~charb() {
    if(_base)
      free(_base);
  }

  size_t capacity() const { return _end - _base; }
  size_t len() const { return _ptr - _base; }

  charb& operator=(const charb &rhs) {
    if(this == &rhs)
      return *this;
    ensure(rhs.capacity());
    memcpy(_base, rhs._base, capacity());
  }

  template<typename T>
  char& operator[](T i) { return _base[i]; }
  template<typename T>
  char operator[](T i) const { return _base[i]; }
  char& operator*() { return *_base; }
  char operator*() const { return *_base; }
  operator char *() const { return _base; }
  operator const char *() const { return _base; }
  operator void *() const { return (void*)_base; }
  bool operator==(const charb &rhs) const { return _base == rhs._base; }
  bool operator!=(const charb &rhs) const { return _base != rhs._base; }
  template<typename T>
  bool operator==(T *rhs) const { return (T*)_base == rhs; }
  template<typename T>
  bool operator!=(T *rhs) const { return (T*)_base != rhs; }
  bool operator!() const { return _base != 0; }

  void ensure(size_t s) {
    size_t clen = _end - _base;
    if(s <= clen)
      return;
    if(s <= 2 * clen)
      s = 2 * clen;
    if(s == 0)
      s = 1024;
    char *nbase = (char *)realloc(_base, s);
    if(!nbase)
      throw std::runtime_error("Error allocating memory");
    _ptr  = nbase + (_ptr - _base);
    _end  = nbase + s;
    _base = nbase;
  }
  void enlarge() { ensure(capacity() * 2); }

  friend char *fgets(charb &b, FILE *stream);
  friend int vsprintf(charb &b, const char *format, va_list ap);
  friend ssize_t getline(charb &b, FILE *stream);
  friend ssize_t getdelim(charb &b, int delim, FILE *stream);
};

/** fgets for char buffer. Expand the size of the buffer if the line
 * does not fit.
 */
char *fgets(charb &b, FILE *stream) {
  char *cptr = b._base;
  long  pos  = ftell(stream);
  long  npos = pos;
  
  while(true) {
    char *res = fgets(cptr, b.capacity() - (cptr - b._base), stream);
    if(!res)
      break;
    size_t char_read;
    if(pos == -1) {
      char_read = strlen(res);
    } else {
      npos = ftell(stream);
      char_read = npos - pos;
      pos = npos;
    }
    cptr      += char_read;
    if(cptr < b._end - 1 || *(cptr - 1) == '\n')
      break;
    size_t off = cptr - b._base;
    b.enlarge();
    cptr = b._base + off;
  }
  
  if(cptr == b._base)
    return 0;
  b._ptr = cptr;
  return b._base;
}

/** Backward compatible fgets for char buffer. The size argument is ignored and present
 * only for backward compatibility.
 */
template<typename T>
char *fgets(charb &b, T size, FILE *stream) { return fgets(b, stream); }

/** Getline for char buffer.
 */
ssize_t getline(charb &b, FILE *stream) {
  size_t n = b.capacity();
  ssize_t res = getline(&b._base, &n, stream);
  if(res == -1)
    return res;
  b._ptr = b._base + res;
  if(n != b.capacity())
    b._end = b._base + n;

  return res;
}

template<typename T>
ssize_t getline(charb &b, T *n, int delim, FILE *stream) { return getline(b, stream); }

/** Getdelim for char buffer.
 */
ssize_t getdelim(charb &b, int delim, FILE *stream) {
  size_t n = b.capacity();
  ssize_t res = getdelim(&b._base, &n, delim, stream);
  if(res == -1)
    return res;
  b._ptr = b._base + res;
  if(n != b.capacity())
    b._end = b._base + n;

  return res;
}

template<typename T>
ssize_t getdelim(charb &b, T *n, int delim, FILE *stream) { return getdelim(b, delim, stream); }


/** Sprintf for buffer. The buffer grows as needed. The return value
 * is < 0 in case of error, or the number of characters written to the
 * char buffer.
 */
int sprintf(charb &b, const char *format, ...) {
  va_list ap;

  va_start(ap, format);
  int res = vsprintf(b, format, ap);
  va_end(ap);

  return res;
}

/** Snprintf for backward compatibility.
 */
template<typename T>
int snprintf(charb &b, T size, const char *format, ...) {
  va_list ap;
  va_start(ap, format);
  int res = vsprintf(b, format, ap);
  va_end(ap);

  return res;
}

/** Vsnprintf for backward compatibility.
 */
template<typename T>
int vsnprintf(charb &b, T size, const char *format, va_list ap) {
  return vsprintf(b, format, ap);
}

int vsprintf(charb &b, const char *format, va_list _ap) {
  int res;
  while(true) {
    va_list ap;
    va_copy(ap, _ap);
    res = vsnprintf(b._base, b.capacity(), format, ap);
    va_end(ap);
    if(res < 0)
      return res;
    if((size_t)res < b.capacity())
      break;
    b.ensure(res + 1);
  }
  b._ptr = b._base + res;
  return res;
}
#endif /* __CHARB_HPP__ */
