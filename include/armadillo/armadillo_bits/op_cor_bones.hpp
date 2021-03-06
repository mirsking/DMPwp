// Copyright (C) 2009-2010 NICTA (www.nicta.com.au)
// Copyright (C) 2009-2010 Conrad Sanderson
// Copyright (C) 2009-2010 Dimitrios Bouzas
// 
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.



//! \addtogroup op_cor
//! @{



class op_cor
  {
  public:
  
  template<typename eT> inline static void direct_cor(Mat<eT>&                out, const Mat<eT>& X,                const uword norm_type);
  template<typename  T> inline static void direct_cor(Mat< std::complex<T> >& out, const Mat< std::complex<T> >& X, const uword norm_type);
  
  template<typename T1> inline static void apply(Mat<typename T1::elem_type>& out, const Op<T1,op_cor>& in);
  };



//! @}
