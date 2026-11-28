/*!
        @file    source.h

        @brief

        @author  Tatsumi Aoyama  (aoym)
                 $LastChangedBy: matufuru $

        @date    $LastChangedDate:: 2023-03-20 10:52:44 #$

        @version $LastChangedRevision: 2499 $
*/

#ifndef SOURCE_INCLUDED
#define SOURCE_INCLUDED

#include "Parameters/parameters.h"
#include "Field/field.h"

#include "IO/bridgeIO.h"

#ifdef USE_FACTORY
#include "Tools/factory.h"
#endif

//! Base class of source for a linear solver.

/*!
    (Coding history will be recovered from trac.)
    Parameters_Source_All is implemented for source_selector.
                                      [ 2 Feb 2013 Y.Namekawa]
    Bug fix: initialization of verbose level was added.
                                      [23 May 2016 H.Matsufuru]
    Add set_all_color_spin,etc        [ 4 Apr 2017 Y.Namekawa]
 */

class Source
{
 public:
  Source() {}
  virtual ~Source() {}

 private:
  // non-copyable
  Source(const Source&);
  Source& operator=(const Source&);

 public:
  virtual void set_parameters(const Parameters&) = 0;

  virtual void get_parameters(Parameters&) const = 0;

  virtual void set(Field&, const int)            = 0;
  virtual void set(Field&, const int, const int) = 0;
  virtual void set_all_color(Field&, const int)  = 0;
  virtual void set_all_color_spin(Field&)        = 0;

#ifdef USE_FACTORY
 public:
  typedef Source *(*ProductCreator)();
  typedef Source *(*ProductCreator_params)(const Parameters&);

  typedef FactoryTemplate<Source, ProductCreator>          Factory;
  typedef FactoryTemplate<Source, ProductCreator_params>   Factory_params;

  static Source *New(const IdentifierType& subtype)
  {
    ProductCreator p = Factory::Find(subtype);
    return p ? (*p)() : 0;
  }

  static Source *New(const IdentifierType& subtype, const Parameters& params)
  {
    ProductCreator_params p = Factory_params::Find(subtype);
    return p ? (*p)(params) : 0;
  }

#ifdef USE_FACTORY_AUTOREGISTER
#else
  static bool init_factory();
#endif
#endif
};
/*!
        @file    source_Local.h

        @brief

        @author  Hideo Matsufuru (matsufuru)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Local source for 4-spinor fermion.

/*!
    This class sets an local source vector for the 4-spinor
    (Wilson-type) fermion.
                                     [19 Feb 2012 H.Matsufuru]
    (Coding history will be recovered from trac.)
    YAML is implemented.             [14 Nov 2012 Y.Namekawa]
    Add set_all_color_spin,etc       [ 4 Apr 2017 Y.Namekawa]
 */

class Source_Local : public Source {
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  Index_lex m_index;          // lexical only.
  std::vector<int> m_source_position;
  bool m_in_node;

 public:
  Source_Local() : m_vl(CommonParameters::Vlevel()) {}

  Source_Local(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    set_parameters(params);
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const std::vector<int>& source_position);

  void get_parameters(Parameters& params) const;

  void set(Field& v, const int idx);
  void set(Field& v, const int i_color, const int i_spin);
  void set_all_color(Field& v, const int i_spin);
  void set_all_color_spin(Field& v);

#ifdef USE_FACTORY
 private:
  static Source *create_object()
  {
    return new Source_Local();
  }

  static Source *create_object_with_params(const Parameters& params)
  {
    return new Source_Local(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Source::Factory::Register("Local", create_object);
    init &= Source::Factory_params::Register("Local", create_object_with_params);
    return init;
  }
#endif
};

/*!
        @file    source_Wall.h

        @brief

        @author  Yusuke Namekawa (namekawa)
                 $LastChangedBy: aoyama $

        @date    $LastChangedDate:: 2021-09-08 06:55:16 #$

        @version $LastChangedRevision: 2314 $
*/

#include "Field/index_lex.h"

#include "IO/bridgeIO.h"
using Bridge::vout;

//! Wall source for 4-spinor fermion.

/*!
    This class sets a wall source vector
    for the 4-spinor (Wilson-type) fermion.
                                [02 Feb 2013 Y.Namekawa]
    Add set_all_color_spin,etc  [ 4 Apr 2017 Y.Namekawa]
 */

class Source_Wall : public Source
{
 public:
  static const std::string class_name;

 private:
  Bridge::VerboseLevel m_vl;

  Index_lex m_index;
  std::vector<int> m_source_position;
  bool m_in_node;

 public:
  Source_Wall() : m_vl(CommonParameters::Vlevel()) {}

  Source_Wall(const Parameters& params)
    : m_vl(CommonParameters::Vlevel())
  {
    set_parameters(params);
  }

  void set_parameters(const Parameters& params);
  void set_parameters(const std::vector<int>& source_position);

  void get_parameters(Parameters& params) const;

  void set(Field& v, const int idx);
  void set(Field& v, const int i_color, const int i_spin);
  void set_all_color(Field& v, const int i_spin);
  void set_all_color_spin(Field& v);

#ifdef USE_FACTORY
 private:
  static Source *create_object()
  {
    return new Source_Wall();
  }

  static Source *create_object_with_params(const Parameters& params)
  {
    return new Source_Wall(params);
  }

 public:
  static bool register_factory()
  {
    bool init = true;
    init &= Source::Factory::Register("Wall", create_object);
    init &= Source::Factory_params::Register("Wall", create_object_with_params);
    return init;
  }
#endif
};

#endif  // SOURCE_INCLUDED
