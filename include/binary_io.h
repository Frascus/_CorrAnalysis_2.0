#ifndef __binary_io__
#define __binary_io__

#include "numerics.h"

  //! binary read, non-vector case
  template <class T>
  auto bin_read(T &out,FILE *file) //const -> enable_if_t<is_pod<T>::value>
  {
    int rc=fread(&out,sizeof(T),1,file);
    if(rc != 1) crash("In bin write, error while reading "+to_string(sizeof(T))+" bytes. Exiting..");
  }

  //! return what read
  template <class T>
  T bin_read()
  {
    T out;
    FILE *file;
    bin_read(out,file);
    return out;
  }

  
  //! binary read, non-vector case
  template <class T>
  auto bin_write(T &out,FILE *file) //const -> enable_if_t<is_pod<T>::value>
  {
    int rc=fwrite(&out,sizeof(T),1,file);
    if(rc != 1) crash("In bin write, error while writing "+to_string(sizeof(T))+" bytes. Exiting..");
  }

  //! return what read
  template <class T>
  T bin_write()
  {
    T out;
    FILE *file;
    bin_write(out,file);
    return out;
  }


















#endif
