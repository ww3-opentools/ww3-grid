/* File: GenCellSidemodule.c
 * This file is auto-generated with f2py (version:2).
 * f2py is a Fortran to Python Interface Generator (FPIG), Second Edition,
 * written by Pearu Peterson <pearu@cens.ioc.ee>.
 * Generation date: Tue Jul 16 04:18:47 2019
 * Do not edit this file directly unless you know what you are doing!!!
 */

#ifdef __cplusplus
extern "C" {
#endif

/*********************** See f2py2e/cfuncs.py: includes ***********************/
#include "Python.h"
#include <stdarg.h>
#include "fortranobject.h"
#include <string.h>

/**************** See f2py2e/rules.py: mod_rules['modulebody'] ****************/
static PyObject *GenCellSide_error;
static PyObject *GenCellSide_module;

/*********************** See f2py2e/cfuncs.py: typedefs ***********************/
typedef char * string;

/****************** See f2py2e/cfuncs.py: typedefs_generated ******************/
/*need_typedefs_generated*/

/********************** See f2py2e/cfuncs.py: cppmacros **********************/
\
#define FAILNULL(p) do {                                            \
    if ((p) == NULL) {                                              \
        PyErr_SetString(PyExc_MemoryError, "NULL pointer found");   \
        goto capi_fail;                                             \
    }                                                               \
} while (0)

#define STRINGMALLOC(str,len)\
    if ((str = (string)malloc(sizeof(char)*(len+1))) == NULL) {\
        PyErr_SetString(PyExc_MemoryError, "out of memory");\
        goto capi_fail;\
    } else {\
        (str)[len] = '\0';\
    }

#define rank(var) var ## _Rank
#define shape(var,dim) var ## _Dims[dim]
#define old_rank(var) (PyArray_NDIM((PyArrayObject *)(capi_ ## var ## _tmp)))
#define old_shape(var,dim) PyArray_DIM(((PyArrayObject *)(capi_ ## var ## _tmp)),dim)
#define fshape(var,dim) shape(var,rank(var)-dim-1)
#define len(var) shape(var,0)
#define flen(var) fshape(var,0)
#define old_size(var) PyArray_SIZE((PyArrayObject *)(capi_ ## var ## _tmp))
/* #define index(i) capi_i ## i */
#define slen(var) capi_ ## var ## _len
#define size(var, ...) f2py_size((PyArrayObject *)(capi_ ## var ## _tmp), ## __VA_ARGS__, -1)

#define STRINGFREE(str) do {if (!(str == NULL)) free(str);} while (0)

#ifdef DEBUGCFUNCS
#define CFUNCSMESS(mess) fprintf(stderr,"debug-capi:"mess);
#define CFUNCSMESSPY(mess,obj) CFUNCSMESS(mess) \
    PyObject_Print((PyObject *)obj,stderr,Py_PRINT_RAW);\
    fprintf(stderr,"\n");
#else
#define CFUNCSMESS(mess)
#define CFUNCSMESSPY(mess,obj)
#endif

#ifndef max
#define max(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef min
#define min(a,b) ((a < b) ? (a) : (b))
#endif
#ifndef MAX
#define MAX(a,b) ((a > b) ? (a) : (b))
#endif
#ifndef MIN
#define MIN(a,b) ((a < b) ? (a) : (b))
#endif

#if defined(PREPEND_FORTRAN)
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F
#else
#define F_FUNC(f,F) _##f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) _##F##_
#else
#define F_FUNC(f,F) _##f##_
#endif
#endif
#else
#if defined(NO_APPEND_FORTRAN)
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F
#else
#define F_FUNC(f,F) f
#endif
#else
#if defined(UPPERCASE_FORTRAN)
#define F_FUNC(f,F) F##_
#else
#define F_FUNC(f,F) f##_
#endif
#endif
#endif
#if defined(UNDERSCORE_G77)
#define F_FUNC_US(f,F) F_FUNC(f##_,F##_)
#else
#define F_FUNC_US(f,F) F_FUNC(f,F)
#endif

#define STRINGCOPYN(to,from,buf_size)                           \
    do {                                                        \
        int _m = (buf_size);                                    \
        char *_to = (to);                                       \
        char *_from = (from);                                   \
        FAILNULL(_to); FAILNULL(_from);                         \
        (void)strncpy(_to, _from, sizeof(char)*_m);             \
        _to[_m-1] = '\0';                                      \
        /* Padding with spaces instead of nulls */              \
        for (_m -= 2; _m >= 0 && _to[_m] == '\0'; _m--) {      \
            _to[_m] = ' ';                                      \
        }                                                       \
    } while (0)


/************************ See f2py2e/cfuncs.py: cfuncs ************************/
static int double_from_pyobj(double* v,PyObject *obj,const char *errmess) {
    PyObject* tmp = NULL;
    if (PyFloat_Check(obj)) {
#ifdef __sgi
        *v = PyFloat_AsDouble(obj);
#else
        *v = PyFloat_AS_DOUBLE(obj);
#endif
        return 1;
    }
    tmp = PyNumber_Float(obj);
    if (tmp) {
#ifdef __sgi
        *v = PyFloat_AsDouble(tmp);
#else
        *v = PyFloat_AS_DOUBLE(tmp);
#endif
        Py_DECREF(tmp);
        return 1;
    }
    if (PyComplex_Check(obj))
        tmp = PyObject_GetAttrString(obj,"real");
    else if (PyString_Check(obj) || PyUnicode_Check(obj))
        /*pass*/;
    else if (PySequence_Check(obj))
        tmp = PySequence_GetItem(obj,0);
    if (tmp) {
        PyErr_Clear();
        if (double_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
        Py_DECREF(tmp);
    }
    {
        PyObject* err = PyErr_Occurred();
        if (err==NULL) err = GenCellSide_error;
        PyErr_SetString(err,errmess);
    }
    return 0;
}

static int f2py_size(PyArrayObject* var, ...)
{
  npy_int sz = 0;
  npy_int dim;
  npy_int rank;
  va_list argp;
  va_start(argp, var);
  dim = va_arg(argp, npy_int);
  if (dim==-1)
    {
      sz = PyArray_SIZE(var);
    }
  else
    {
      rank = PyArray_NDIM(var);
      if (dim>=1 && dim<=rank)
        sz = PyArray_DIM(var, dim-1);
      else
        fprintf(stderr, "f2py_size: 2nd argument value=%d fails to satisfy 1<=value<=%d. Result will be 0.\n", dim, rank);
    }
  va_end(argp);
  return sz;
}

static int string_from_pyobj(string *str,int *len,const string inistr,PyObject *obj,const char *errmess) {
    PyArrayObject *arr = NULL;
    PyObject *tmp = NULL;
#ifdef DEBUGCFUNCS
fprintf(stderr,"string_from_pyobj(str='%s',len=%d,inistr='%s',obj=%p)\n",(char*)str,*len,(char *)inistr,obj);
#endif
    if (obj == Py_None) {
        if (*len == -1)
            *len = strlen(inistr); /* Will this cause problems? */
        STRINGMALLOC(*str,*len);
        STRINGCOPYN(*str,inistr,*len+1);
        return 1;
    }
    if (PyArray_Check(obj)) {
        if ((arr = (PyArrayObject *)obj) == NULL)
            goto capi_fail;
        if (!ISCONTIGUOUS(arr)) {
            PyErr_SetString(PyExc_ValueError,"array object is non-contiguous.");
            goto capi_fail;
        }
        if (*len == -1)
            *len = (PyArray_ITEMSIZE(arr))*PyArray_SIZE(arr);
        STRINGMALLOC(*str,*len);
        STRINGCOPYN(*str,PyArray_DATA(arr),*len+1);
        return 1;
    }
    if (PyString_Check(obj)) {
        tmp = obj;
        Py_INCREF(tmp);
    }
#if PY_VERSION_HEX >= 0x03000000
    else if (PyUnicode_Check(obj)) {
        tmp = PyUnicode_AsASCIIString(obj);
    }
    else {
        PyObject *tmp2;
        tmp2 = PyObject_Str(obj);
        if (tmp2) {
            tmp = PyUnicode_AsASCIIString(tmp2);
            Py_DECREF(tmp2);
        }
        else {
            tmp = NULL;
        }
    }
#else
    else {
        tmp = PyObject_Str(obj);
    }
#endif
    if (tmp == NULL) goto capi_fail;
    if (*len == -1)
        *len = PyString_GET_SIZE(tmp);
    STRINGMALLOC(*str,*len);
    STRINGCOPYN(*str,PyString_AS_STRING(tmp),*len+1);
    Py_DECREF(tmp);
    return 1;
capi_fail:
    Py_XDECREF(tmp);
    {
        PyObject* err = PyErr_Occurred();
        if (err==NULL) err = GenCellSide_error;
        PyErr_SetString(err,errmess);
    }
    return 0;
}

static int int_from_pyobj(int* v,PyObject *obj,const char *errmess) {
    PyObject* tmp = NULL;
    if (PyInt_Check(obj)) {
        *v = (int)PyInt_AS_LONG(obj);
        return 1;
    }
    tmp = PyNumber_Int(obj);
    if (tmp) {
        *v = PyInt_AS_LONG(tmp);
        Py_DECREF(tmp);
        return 1;
    }
    if (PyComplex_Check(obj))
        tmp = PyObject_GetAttrString(obj,"real");
    else if (PyString_Check(obj) || PyUnicode_Check(obj))
        /*pass*/;
    else if (PySequence_Check(obj))
        tmp = PySequence_GetItem(obj,0);
    if (tmp) {
        PyErr_Clear();
        if (int_from_pyobj(v,tmp,errmess)) {Py_DECREF(tmp); return 1;}
        Py_DECREF(tmp);
    }
    {
        PyObject* err = PyErr_Occurred();
        if (err==NULL) err = GenCellSide_error;
        PyErr_SetString(err,errmess);
    }
    return 0;
}

static int float_from_pyobj(float* v,PyObject *obj,const char *errmess) {
    double d=0.0;
    if (double_from_pyobj(&d,obj,errmess)) {
        *v = (float)d;
        return 1;
    }
    return 0;
}


/********************* See f2py2e/cfuncs.py: userincludes *********************/
/*need_userincludes*/

/********************* See f2py2e/capi_rules.py: usercode *********************/


/* See f2py2e/rules.py */
/*eof externroutines*/

/******************** See f2py2e/capi_rules.py: usercode1 ********************/


/******************* See f2py2e/cb_rules.py: buildcallback *******************/
/*need_callbacks*/

/*********************** See f2py2e/rules.py: buildapi ***********************/

/********************************** adapgrid **********************************/
static char doc_f2py_rout_GenCellSide_gencellside_adapgrid[] = "\
adapgrid(gridid,nlat,nlon,bx,by)\n\nWrapper for ``adapgrid``.\
\n\nParameters\n----------\n"
"gridid : input string(len=32)\n"
"nlat : input int\n"
"nlon : input int\n"
"bx : input float\n"
"by : input float";
/*  */
static PyObject *f2py_rout_GenCellSide_gencellside_adapgrid(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(string,int*,int*,float*,float*,size_t)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  string gridid = NULL;
  int slen(gridid);
  PyObject *gridid_capi = Py_None;
  int nlat = 0;
  PyObject *nlat_capi = Py_None;
  int nlon = 0;
  PyObject *nlon_capi = Py_None;
  float bx = 0;
  PyObject *bx_capi = Py_None;
  float by = 0;
  PyObject *by_capi = Py_None;
  static char *capi_kwlist[] = {"gridid","nlat","nlon","bx","by",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OOOOO:GenCellSide.gencellside.adapgrid",\
    capi_kwlist,&gridid_capi,&nlat_capi,&nlon_capi,&bx_capi,&by_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable gridid */
  slen(gridid) = 32;
  f2py_success = string_from_pyobj(&gridid,&slen(gridid),"",gridid_capi,"string_from_pyobj failed in converting 1st argument `gridid' of GenCellSide.gencellside.adapgrid to C string");
  if (f2py_success) {
  /* Processing variable nlat */
    f2py_success = int_from_pyobj(&nlat,nlat_capi,"GenCellSide.gencellside.adapgrid() 2nd argument (nlat) can't be converted to int");
  if (f2py_success) {
  /* Processing variable nlon */
    f2py_success = int_from_pyobj(&nlon,nlon_capi,"GenCellSide.gencellside.adapgrid() 3rd argument (nlon) can't be converted to int");
  if (f2py_success) {
  /* Processing variable bx */
    f2py_success = float_from_pyobj(&bx,bx_capi,"GenCellSide.gencellside.adapgrid() 4th argument (bx) can't be converted to float");
  if (f2py_success) {
  /* Processing variable by */
    f2py_success = float_from_pyobj(&by,by_capi,"GenCellSide.gencellside.adapgrid() 5th argument (by) can't be converted to float");
  if (f2py_success) {
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
        (*f2py_func)(gridid,&nlat,&nlon,&bx,&by,slen(gridid));
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("");
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  } /*if (f2py_success) of by*/
  /* End of cleaning variable by */
  } /*if (f2py_success) of bx*/
  /* End of cleaning variable bx */
  } /*if (f2py_success) of nlon*/
  /* End of cleaning variable nlon */
  } /*if (f2py_success) of nlat*/
  /* End of cleaning variable nlat */
    STRINGFREE(gridid);
  }  /*if (f2py_success) of gridid*/
  /* End of cleaning variable gridid */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/****************************** end of adapgrid ******************************/

/********************************** cellside **********************************/
static char doc_f2py_rout_GenCellSide_gencellside_cellside[] = "\
cellside(gridid,nlon)\n\nWrapper for ``cellside``.\
\n\nParameters\n----------\n"
"gridid : input string(len=32)\n"
"nlon : input int";
/*  */
static PyObject *f2py_rout_GenCellSide_gencellside_cellside(const PyObject *capi_self,
                           PyObject *capi_args,
                           PyObject *capi_keywds,
                           void (*f2py_func)(string,int*,size_t)) {
  PyObject * volatile capi_buildvalue = NULL;
  volatile int f2py_success = 1;
/*decl*/

  string gridid = NULL;
  int slen(gridid);
  PyObject *gridid_capi = Py_None;
  int nlon = 0;
  PyObject *nlon_capi = Py_None;
  static char *capi_kwlist[] = {"gridid","nlon",NULL};

/*routdebugenter*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_clock();
#endif
  if (!PyArg_ParseTupleAndKeywords(capi_args,capi_keywds,\
    "OO:GenCellSide.gencellside.cellside",\
    capi_kwlist,&gridid_capi,&nlon_capi))
    return NULL;
/*frompyobj*/
  /* Processing variable gridid */
  slen(gridid) = 32;
  f2py_success = string_from_pyobj(&gridid,&slen(gridid),"",gridid_capi,"string_from_pyobj failed in converting 1st argument `gridid' of GenCellSide.gencellside.cellside to C string");
  if (f2py_success) {
  /* Processing variable nlon */
    f2py_success = int_from_pyobj(&nlon,nlon_capi,"GenCellSide.gencellside.cellside() 2nd argument (nlon) can't be converted to int");
  if (f2py_success) {
/*end of frompyobj*/
#ifdef F2PY_REPORT_ATEXIT
f2py_start_call_clock();
#endif
/*callfortranroutine*/
        (*f2py_func)(gridid,&nlon,slen(gridid));
if (PyErr_Occurred())
  f2py_success = 0;
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_call_clock();
#endif
/*end of callfortranroutine*/
    if (f2py_success) {
/*pyobjfrom*/
/*end of pyobjfrom*/
    CFUNCSMESS("Building return value.\n");
    capi_buildvalue = Py_BuildValue("");
/*closepyobjfrom*/
/*end of closepyobjfrom*/
    } /*if (f2py_success) after callfortranroutine*/
/*cleanupfrompyobj*/
  } /*if (f2py_success) of nlon*/
  /* End of cleaning variable nlon */
    STRINGFREE(gridid);
  }  /*if (f2py_success) of gridid*/
  /* End of cleaning variable gridid */
/*end of cleanupfrompyobj*/
  if (capi_buildvalue == NULL) {
/*routdebugfailure*/
  } else {
/*routdebugleave*/
  }
  CFUNCSMESS("Freeing memory.\n");
/*freemem*/
#ifdef F2PY_REPORT_ATEXIT
f2py_stop_clock();
#endif
  return capi_buildvalue;
}
/****************************** end of cellside ******************************/
/*eof body*/

/******************* See f2py2e/f90mod_rules.py: buildhooks *******************/

static FortranDataDef f2py_gencellside_def[] = {
  {"ncl",0,{{-1}},NPY_INT},
  {"nfc",0,{{-1}},NPY_INT},
  {"dx",1,{{2250001}},NPY_FLOAT},
  {"dy",1,{{2250001}},NPY_FLOAT},
  {"dxr",1,{{2250001}},NPY_FLOAT},
  {"dyr",1,{{2250001}},NPY_FLOAT},
  {"chg",1,{{2250001}},NPY_FLOAT},
  {"nu",0,{{-1}},NPY_INT},
  {"nv",0,{{-1}},NPY_INT},
  {"nc",0,{{-1}},NPY_INT},
  {"ns",0,{{-1}},NPY_INT},
  {"nt",0,{{-1}},NPY_INT},
  {"n9",0,{{-1}},NPY_INT},
  {"n8",0,{{-1}},NPY_INT},
  {"n4",0,{{-1}},NPY_INT},
  {"n2",0,{{-1}},NPY_INT},
  {"n1",0,{{-1}},NPY_INT},
  {"ice",2,{{4,2250000}},NPY_INT},
  {"kg",1,{{2250000}},NPY_INT},
  {"isd",2,{{7,2360000}},NPY_INT},
  {"jsd",2,{{8,2360000}},NPY_INT},
  {"i",0,{{-1}},NPY_INT},
  {"ii",0,{{-1}},NPY_INT},
  {"ij",0,{{-1}},NPY_INT},
  {"ijk",0,{{-1}},NPY_INT},
  {"j",0,{{-1}},NPY_INT},
  {"jj",0,{{-1}},NPY_INT},
  {"jk",0,{{-1}},NPY_INT},
  {"jkl",0,{{-1}},NPY_INT},
  {"k",0,{{-1}},NPY_INT},
  {"kk",0,{{-1}},NPY_INT},
  {"kl",0,{{-1}},NPY_INT},
  {"klm",0,{{-1}},NPY_INT},
  {"l",0,{{-1}},NPY_INT},
  {"ll",0,{{-1}},NPY_INT},
  {"lm",0,{{-1}},NPY_INT},
  {"lmn",0,{{-1}},NPY_INT},
  {"m",0,{{-1}},NPY_INT},
  {"mm",0,{{-1}},NPY_INT},
  {"mn",0,{{-1}},NPY_INT},
  {"n",0,{{-1}},NPY_INT},
  {"nn",0,{{-1}},NPY_INT},
  {"nlat",0,{{-1}},NPY_INT},
  {"nlon",0,{{-1}},NPY_INT},
  {"xext",1,{{6,1}},NPY_STRING},
  {"rundate",1,{{16}},NPY_STRING},
  {"adapgrid",-1,{{-1}},0,NULL,(void *)f2py_rout_GenCellSide_gencellside_adapgrid,doc_f2py_rout_GenCellSide_gencellside_adapgrid},
  {"cellside",-1,{{-1}},0,NULL,(void *)f2py_rout_GenCellSide_gencellside_cellside,doc_f2py_rout_GenCellSide_gencellside_cellside},
  {NULL}
};

static void f2py_setup_gencellside(char *ncl,char *nfc,char *dx,char *dy,char *dxr,char *dyr,char *chg,char *nu,char *nv,char *nc,char *ns,char *nt,char *n9,char *n8,char *n4,char *n2,char *n1,char *ice,char *kg,char *isd,char *jsd,char *i,char *ii,char *ij,char *ijk,char *j,char *jj,char *jk,char *jkl,char *k,char *kk,char *kl,char *klm,char *l,char *ll,char *lm,char *lmn,char *m,char *mm,char *mn,char *n,char *nn,char *nlat,char *nlon,char *xext,char *rundate,char *adapgrid,char *cellside) {
  int i_f2py=0;
  f2py_gencellside_def[i_f2py++].data = ncl;
  f2py_gencellside_def[i_f2py++].data = nfc;
  f2py_gencellside_def[i_f2py++].data = dx;
  f2py_gencellside_def[i_f2py++].data = dy;
  f2py_gencellside_def[i_f2py++].data = dxr;
  f2py_gencellside_def[i_f2py++].data = dyr;
  f2py_gencellside_def[i_f2py++].data = chg;
  f2py_gencellside_def[i_f2py++].data = nu;
  f2py_gencellside_def[i_f2py++].data = nv;
  f2py_gencellside_def[i_f2py++].data = nc;
  f2py_gencellside_def[i_f2py++].data = ns;
  f2py_gencellside_def[i_f2py++].data = nt;
  f2py_gencellside_def[i_f2py++].data = n9;
  f2py_gencellside_def[i_f2py++].data = n8;
  f2py_gencellside_def[i_f2py++].data = n4;
  f2py_gencellside_def[i_f2py++].data = n2;
  f2py_gencellside_def[i_f2py++].data = n1;
  f2py_gencellside_def[i_f2py++].data = ice;
  f2py_gencellside_def[i_f2py++].data = kg;
  f2py_gencellside_def[i_f2py++].data = isd;
  f2py_gencellside_def[i_f2py++].data = jsd;
  f2py_gencellside_def[i_f2py++].data = i;
  f2py_gencellside_def[i_f2py++].data = ii;
  f2py_gencellside_def[i_f2py++].data = ij;
  f2py_gencellside_def[i_f2py++].data = ijk;
  f2py_gencellside_def[i_f2py++].data = j;
  f2py_gencellside_def[i_f2py++].data = jj;
  f2py_gencellside_def[i_f2py++].data = jk;
  f2py_gencellside_def[i_f2py++].data = jkl;
  f2py_gencellside_def[i_f2py++].data = k;
  f2py_gencellside_def[i_f2py++].data = kk;
  f2py_gencellside_def[i_f2py++].data = kl;
  f2py_gencellside_def[i_f2py++].data = klm;
  f2py_gencellside_def[i_f2py++].data = l;
  f2py_gencellside_def[i_f2py++].data = ll;
  f2py_gencellside_def[i_f2py++].data = lm;
  f2py_gencellside_def[i_f2py++].data = lmn;
  f2py_gencellside_def[i_f2py++].data = m;
  f2py_gencellside_def[i_f2py++].data = mm;
  f2py_gencellside_def[i_f2py++].data = mn;
  f2py_gencellside_def[i_f2py++].data = n;
  f2py_gencellside_def[i_f2py++].data = nn;
  f2py_gencellside_def[i_f2py++].data = nlat;
  f2py_gencellside_def[i_f2py++].data = nlon;
  f2py_gencellside_def[i_f2py++].data = xext;
  f2py_gencellside_def[i_f2py++].data = rundate;
  f2py_gencellside_def[i_f2py++].data = adapgrid;
  f2py_gencellside_def[i_f2py++].data = cellside;
}
extern void F_FUNC(f2pyinitgencellside,F2PYINITGENCELLSIDE)(void (*)(char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char*,char *,char *));
static void f2py_init_gencellside(void) {
  F_FUNC(f2pyinitgencellside,F2PYINITGENCELLSIDE)(f2py_setup_gencellside);
}

/*need_f90modhooks*/

/************** See f2py2e/rules.py: module_rules['modulebody'] **************/

/******************* See f2py2e/common_rules.py: buildhooks *******************/

/*need_commonhooks*/

/**************************** See f2py2e/rules.py ****************************/

static FortranDataDef f2py_routine_defs[] = {

/*eof routine_defs*/
  {NULL}
};

static PyMethodDef f2py_module_methods[] = {

  {NULL,NULL}
};

#if PY_VERSION_HEX >= 0x03000000
static struct PyModuleDef moduledef = {
  PyModuleDef_HEAD_INIT,
  "GenCellSide",
  NULL,
  -1,
  f2py_module_methods,
  NULL,
  NULL,
  NULL,
  NULL
};
#endif

#if PY_VERSION_HEX >= 0x03000000
#define RETVAL m
PyMODINIT_FUNC PyInit_GenCellSide(void) {
#else
#define RETVAL
PyMODINIT_FUNC initGenCellSide(void) {
#endif
  int i;
  PyObject *m,*d, *s;
#if PY_VERSION_HEX >= 0x03000000
  m = GenCellSide_module = PyModule_Create(&moduledef);
#else
  m = GenCellSide_module = Py_InitModule("GenCellSide", f2py_module_methods);
#endif
  Py_TYPE(&PyFortran_Type) = &PyType_Type;
  import_array();
  if (PyErr_Occurred())
    {PyErr_SetString(PyExc_ImportError, "can't initialize module GenCellSide (failed to import numpy)"); return RETVAL;}
  d = PyModule_GetDict(m);
  s = PyString_FromString("$Revision: $");
  PyDict_SetItemString(d, "__version__", s);
#if PY_VERSION_HEX >= 0x03000000
  s = PyUnicode_FromString(
#else
  s = PyString_FromString(
#endif
    "This module 'GenCellSide' is auto-generated with f2py (version:2).\nFunctions:\n"
"Fortran 90/95 modules:\n""  gencellside --- ncl,nfc,dx,dy,dxr,dyr,chg,nu,nv,nc,ns,nt,n9,n8,n4,n2,n1,ice,kg,isd,jsd,i,ii,ij,ijk,j,jj,jk,jkl,k,kk,kl,klm,l,ll,lm,lmn,m,mm,mn,n,nn,nlat,nlon,xext,rundate,adapgrid(),cellside()"".");
  PyDict_SetItemString(d, "__doc__", s);
  GenCellSide_error = PyErr_NewException ("GenCellSide.error", NULL, NULL);
  Py_DECREF(s);
  for(i=0;f2py_routine_defs[i].name!=NULL;i++)
    PyDict_SetItemString(d, f2py_routine_defs[i].name,PyFortranObject_NewAsAttr(&f2py_routine_defs[i]));


/*eof initf2pywraphooks*/
  PyDict_SetItemString(d, "gencellside", PyFortranObject_New(f2py_gencellside_def,f2py_init_gencellside));
/*eof initf90modhooks*/

/*eof initcommonhooks*/


#ifdef F2PY_REPORT_ATEXIT
  if (! PyErr_Occurred())
    on_exit(f2py_report_on_exit,(void*)"GenCellSide");
#endif

  return RETVAL;
}
#ifdef __cplusplus
}
#endif
