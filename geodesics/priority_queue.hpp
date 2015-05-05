/** 
 * \file priority_queue.hpp
 *  
 * \brief    Definition    and    implementation    of    the    class
 * dynamic_priority_queue,  which allows for  element deletion  in \f$
 * O(n \log n)\f$ time, where \f$n\f$ is the number of elements in the
 * queue.
 *
 * \author
 * Marcelo Ferreira Siqueira \n
 * Universidade Federal do Rio Grande do Norte, \n
 * Departamento de Informatica e Matematica Aplicada, \n
 * marcelo at dimap (dot) ufrn (dot) br
 *
 * \version 1.0
 * \date January 2010
 *
 * \attention This program is distributed WITHOUT ANY WARRANTY, and it
 *            may be freely redistributed under the condition that the
 *            copyright notices  are not removed,  and no compensation
 *            is received. Private, research, and institutional use is
 *            free. Distribution of this  code as part of a commercial
 *            system  is permissible ONLY  BY DIRECT  ARRANGEMENT WITH
 *            THE AUTHOR.
 */

#ifndef PRIORITY_QUEUE_HPP
#define PRIORITY_QUEUE_HPP

#include <algorithm>     // std::make_heap
#include <vector>        // std::vector
#include <functional>    // std::less 
#include <map>           // std::map
#include <cassert>       // assert



/** 
 * \defgroup COMMONNameSpace Namespace common.
 * @{
 */

/**
 * \namespace common
 *
 * \brief   The   namespace  common   contains   the  definition   and
 * implementation  of several  general classes,  which  contains basic
 * algorithms and/or  basic data structures  that can be used  by more
 * specialized classes.
 */

namespace common {

  /**
   * \struct KeyCompare
   *
   * \brief A template-based struct to compare two generic objects.
   */
  template < class Node, class key_compare >
  struct KeyCompare
  {
    // ---------------------------------------------------------------
    //
    // Public methods 
    //
    // ---------------------------------------------------------------
  
    /**
     * \fn bool operator()( const Node& n1 , const Node& n2 )
     *
     * \brief Compares two generic objects and returns the logic value
     * true if they are the same and the logic value false otherwise.
     *
     * \param n1 A reference to a generic object.
     * \param n2 A reference to another generic object.
     *
     * \returns A logic value true if the two given objects are equal,
     * and false otherwise.
     */
    bool operator()( const Node& n1 , const Node& n2 )
    {
      key_compare comp;
      return comp( n1.second , n2.second ) ;
    }

  } ;


  /**
   * \class dynamic_priority_queue 
   *
   * \brief A template-based class  that represents a generic, dynamic
   * priority queue.
   */
  template < 
             typename object_type, 
             typename key_type = object_type, 
             class key_compare = std::less< key_type >
           >
  class  dynamic_priority_queue {
  private:
    // ---------------------------------------------------------------
    //
    // Type name definitions.
    //
    // ---------------------------------------------------------------

    /**
     * \typedef node
     *
     * \brief  Definition  of  a  type  name for  the  priority  queue
     * elements.
     */
    typedef std::pair< object_type , key_type > node ;

    /**
     * \typedef idMap
     *
     * \brief Definition of a type name for a hash table to store
     * indices.
     */
    typedef std::map< object_type , size_t >  idMap ;

    /**
     * \typedef idIterator
     *
     * \brief Definition of a type name for a hash table to store
     * iterators.
     */
    typedef typename std::map< object_type , size_t >::iterator idIterator ;
         
  public:
    // ---------------------------------------------------------------
    //
    // Public methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn dynamic_priority_queue()
     *
     * \brief Creates an instance of this class (an empty heap).
     *
     */    
    dynamic_priority_queue() ;

    /**
     * \fn  dynamic_priority_queue( std::vector< node >& v )
     *
     * \brief Creates an instance of this class.
     *
     * \param v A reference to a vector of elements used to initialize
     * the priority queue.
     */    
    dynamic_priority_queue( std::vector< node >& v ) ;


    /**
     * \fn ~dynamic_priority_queue()
     *
     * \brief Destroys an instance of this class.
     */    
    virtual ~dynamic_priority_queue()
    {}


    /**
     * \fn void change_key_at( size_t at , key_type k )
     *
     * \brief Change the key of  the priority queue element in a given
     * position.
     *
     * \param at Position of the element in the priority queue.
     * \param k Value to be assigned with the element key.
     */    
    void change_key_at( size_t at , key_type k ) ;


    /**
     * \fn void change_key_at( object_type at , key_type k )
     *
     * \brief Change the key of a given priority queue element. If the
     * element  is not  in the  queue, the  program is  aborted  by an
     * exception.
     *
     * \param at Element in the priority queue.
     * \param k Value to be assigned to the element key.
     */    
    void change_key_at( object_type at , key_type k ) ;


    /**
     * \fn inline const object_type& top() const
     *
     * \brief Gets the element on the top of the priority queue.
     *
     * \return A constant  reference to the element on  the top of the
     * priority queue.
     */
    inline const object_type& top() const 
    {
      //
      // The heap cannot be empty.
      //
      assert( !empty() ) ;

      return _h[ 0 ].first ; 
    }


    /**
     * \fn inline const key_type& top_key() const
     *
     * \brief Gets the  key of the element on the  top of the priority
     * queue.
     *
     * \return A constant  reference to the key of  the element on the
     * top of the priority queue.
     */
    inline const key_type& top_key() const 
    { 
      //
      // The heap cannot be empty.
      //

      return _h[ 0 ].second ; 
    }


    /**
     * \fn void pop()
     *
     * \brief Remove the element with highest priority from the heap.
     */
    void pop() ;


    /**
     * \fn void push( object_type obj , key_type key ) 
     *
     * \brief Insert an element into the priority queue.
     *
     * \param obj An  element to be inserted into  the priority queue.
     * \param key  The key  assigned with the  element to  be inserted
     * into the priority queue.
     */    
    void push( object_type , key_type ) ;


    /**
     * \fn void erase( size_t at )
     *
     * \brief Remove an element from the priority queue.
     *
     * \param  at The  index  of an  element  to be  removed from  the
     * priority queue.
     */
    void erase( size_t at ) ; 


    /**
     * \fn void erase( object_type obj )
     *
     * \brief Remove an element from the priority queue.
     *
     * \param obj An element to be removed from the priority queue.
     */
    void erase( object_type at ) ;


    /**
     * \fn inline void clear()
     *
     * \brief Remove all elements from the priority queue.
     */
    inline void clear() 
    {
      _hsize = 0 ;
      _id.clear() ; 
    }
         
     
    /**
     * \fn inline bool empty() const
     *
     * \brief Decides whether the heap is empty.
     *
     * \return  The longic value  true if  the heap  is empty,  and the
     * logic value false otherwise.
     */
    inline bool empty() const 
    { 
      return ( _hsize == 0 ) ; 
    }


    /**
     * \fn inline sizet size() const
     *
     * \brief Computes the number of elements of the priority queue.
     *
     * \return The number of elements of the priority queue.
     */
    size_t size() const 
    { 
      return _hsize ; 
    }


    /**
     * \fn inline bool has(object_type o)
     *
     * \brief Search for an element in the priority queue.
     *
     * \return The longic  value true if the element  is found and the
     * logic value false otherwise.
     */
    inline bool has( object_type obj )
    {
      idIterator it = _id.find( obj ) ;

      return ( it != _id.end() ) ;
    }

  private:
    // ---------------------------------------------------------------
    //
    // Private methods
    //
    // ---------------------------------------------------------------

    /**
     * \fn void go_up( size_t idx )
     *
     * \brief Try to move a given element up in the priority queue.
     *
     * \param idx The index of an element in the priority queue.
     */
    void go_up( size_t idx ) ;


    /**
     * \fn void go_down( size_t idx )
     *
     * \brief Try to move a given element down in the priority queue.
     *
     * \param idx The index of an element in the priority queue.
     */
    void go_down( size_t idx ) ;


    // ---------------------------------------------------------------
    //
    // Private data members.
    //
    // ---------------------------------------------------------------
    std::vector< node > _h ; ///< The heap container
    size_t _hsize ;  ///< Current heap size
    idMap _id ;  ///< Index map
    KeyCompare< node , key_compare > _compare ;  ///< A comparator element

  } ;


  /**
   * \fn  dynamic_priority_queue<object_type, key_type, key_compare>::dynamic_priority_queue()
   *
   * \brief Creates an instance of this class (an empty heap).
   *
   */    
   template < 
             typename object_type , 
             typename key_type ,
             class key_compare
            >
   dynamic_priority_queue<object_type, key_type, key_compare>::dynamic_priority_queue()
   {
     _hsize = 0 ;
     _h.reserve( 100 ) ;

      return ;
   }


  /**
   * \fn  dynamic_priority_queue<object_type, key_type, key_compare>::dynamic_priority_queue( std::vector< node >& v )
   *
   * \brief Creates an instance of this class.
   *
   * \param v A  reference to a vector of  elements used to initialize
   * the priority queue.
   */    
   template < 
             typename object_type , 
             typename key_type ,
             class key_compare
            >
   dynamic_priority_queue<object_type, key_type, key_compare>::dynamic_priority_queue( 
     std::vector< node >& v )
   :
     _h( v ) , _hsize( _h.size() )
   {
      std::make_heap( _h.begin() , _h.end() , _compare ) ;

      for ( size_t i = 0 ; i < _h.size(); i++ ) {
	_id[ _h[ i ].first ] = i ;
      }

      return ;
   }


  /**
   * \fn void dynamic_priority_queue<object_type, key_type, key_compare>::change_key_at( size_t at , key_type k )
   *
   * \brief Change  the key of the  priority queue element  in a given
   * position.
   *
   * \param at Position of the element in the priority queue.
   * \param k Value to be assigned with the element key.
   */
  template < 
             typename object_type , 
             typename key_type ,
             class key_compare
            >
  void
  dynamic_priority_queue< object_type , key_type, key_compare >::change_key_at(
									       size_t at , 
									       key_type k
									      )
  {
    /*
     * Is the index valid?
     */
    assert ( at < _hsize ) ;

    /*
     * If the key value does not change, do nothing. Otherwise, update
     * the priority queue.
     */
    if ( _h[ at ].second != k ) {
      if ( _compare( _h[ at ] , std::make_pair( _h[ at ].first , k ) ) ) {
	_h[ at ].second = k ;
	go_up( at ) ;
      }
      else {
	_h[ at ].second = k ;
        go_down( at ) ;
      }
    }

    return ;
  }


  /**
   * \fn void dynamic_priority_queue<object_type, key_type, key_compare>::change_key_at( object_type at , key_type k )
   *
   * \brief Change the  key of a given priority  queue element. If the
   * element  is  not in  the  queue, the  program  is  aborted by  an
   * exception.
   *
   * \param at Element in the priority queue.
   * \param k Value to be assigned to the element key.
   */    
  template < 
            typename object_type ,
            typename key_type ,
            class key_compare
          >
  void 
  dynamic_priority_queue< object_type, key_type , key_compare >::change_key_at(
									       object_type at ,
									       key_type k
									      )
  {
    idIterator it = _id.find( at ) ;

    assert( it != _id.end() ) ;

    change_key_at( it->second , k ) ;
 
    return ;
  }


  /**
   * \fn void dynamic_priority_queue<object_type, key_type, key_compare>::pop()
   *
   * \brief Remove the element with highest priority from the heap.
   */
   template <
             typename object_type ,
             typename key_type ,
             class key_compare
            >
   void 
   dynamic_priority_queue< object_type , key_type , key_compare >::pop()
   {
     /*
      * Heap cannot be empty!
      */
     assert( !empty() ) ;

     /* Drop first element from id map. */
     _id.erase( _h[ 0 ].first ) ;

     /* Copy the element to the top most position. */
     --_hsize ;
     _h[ 0 ] = _h[ _hsize ] ;

     /* Update id  */
     _id[ _h[ 0 ].first ] = 0 ;

     /* Move down the element in the top most position (if needed). */
     go_down( 0 ) ;

     return ;
   }


   /**
    * \fn void dynamic_priority_queue<object_type, key_type, key_compare>::push( object_type obj , key_type key ) 
    *
    * \brief Insert an element into the priority queue.
    *
    * \param obj An  element to be inserted into  the priority queue.
    * \param key  The key  assigned with the  element to  be inserted
    * into the priority queue.
    */
   template<
            typename object_type ,
            typename key_type ,
            class key_compare
           >
   void 
   dynamic_priority_queue< object_type , key_type , key_compare >::push(
     object_type obj ,
     key_type key 
   )
   {
     if ( _hsize != _h.size() ) {
       _h[ _hsize ] = std::make_pair( obj , key ) ;
     }
     else {
       _h.push_back( std::make_pair( obj , key ) ) ;
     }

     _id[ obj ] = _hsize ;

     go_up( _hsize ) ;

     ++_hsize ;

     return ;
   }


  /**
   * \fn void dynamic_priority_queue<object_type, key_type, key_compare>::erase( size_t at )
   *
   * \brief Remove an element from the priority queue.
   *
   * \param at The index of an element to be removed from the priority
   * queue.
   */
  template < 
            typename object_type ,
            typename key_type ,
            class key_compare
           >
  void
  dynamic_priority_queue< object_type , key_type , key_compare >::erase( size_t at )
  {
    /*
     * Is the given index valid?
     */
    assert( at < _hsize ) ;

    /*
     * Drop first element from id map.
     */
    _id.erase( _h[ at ].first ) ;

    --_hsize ;

    if ( at == _hsize ) {
      return ;
    }

    /*
     * Overwrite the element with index "at" with the last element.
     */
    key_type key1 = _h[     at].second ;
    key_type key2 = _h[ _hsize].second ;

    _h[ at ] = _h[ _hsize ] ;
    _h[ at ].second = key1 ;

    _id[ _h[ at ].first ] = at;

    change_key_at( at , key2 ) ;

    return ;
  }


  /**
   * \fn void dynamic_priority_queue<object_type, key_type, key_compare>::erase( object_type obj )
   *
   * \brief Remove an element from the priority queue.
   *
   * \param obj An element to be removed from the priority queue.
   */
  template <
            typename object_type ,
            typename key_type ,
            class key_compare
           >
  void
  dynamic_priority_queue< object_type , key_type , key_compare >::erase( object_type obj )
  {
    idIterator it = _id.find( obj ) ;

    /*
     * Element must belong to the heap.
     */
    assert( it != _id.end() ) ;

    erase( it->second ) ;

    return ;
  }


  /**
   * \fn void dynamic_priority_queue<object_type, key_type, key_compare>::go_up( size_t idx )
   *
   * \brief Try to move a given element up in the priority queue.
   *
   * \param idx The index of an element in the priority queue.
   */
  template < typename object_type , typename key_type , class key_compare >
  void
  dynamic_priority_queue< object_type , key_type , key_compare >::go_up( size_t idx )
  {
    size_t i = idx ;

    while ( i > 0 ) {
      size_t p = ( i - 1 ) >> 1 ;

      if ( _compare( _h[ p ] , _h[ i ] ) ) {
	node temp = _h[ i ] ;
	_h[ i ] = _h[ p ] ;
	_id[ _h[ i ].first ] = i ;
	_h[ p ] = temp ;
	_id[ _h[ p ].first ] = p ;
      }
      else {
	break ;
      }

      i = p ;
    }

    return ;
  }


  /**
   * \fn void dynamic_priority_queue<object_type, key_type, key_compare>::go_down( size_t idx )
   *
   * \brief Try to move a given element down in the priority queue.
   *
   * \param idx The index of an element in the priority queue.
   */
  template <
            typename object_type ,
            typename key_type ,
            class key_compare
           >
  void
  dynamic_priority_queue< object_type , key_type , key_compare >::go_down( size_t idx )
  {
    size_t i = idx ;

    while ( i < _hsize ) {
      size_t l = 2 * i + 1 ;
      if ( l < _hsize ) {
	if ( ( l + 1 ) < _hsize ) {
	  if ( _compare( _h[ l ], _h[ l + 1 ] ) ) {
	    ++l ;
	  }
	}
	if ( _compare( _h[ i ] , _h[ l ] ) ) {
	  node temp = _h[ i ] ;
	  _h[ i ] = _h[ l ] ;
	  _id[ _h[ i ].first ] = i ;
	  _h[ l ] = temp ;
	  _id[ _h[ l ].first] = l ;

	  i = l ;
	}
	else {
	  i = _hsize ;
	}
      }
      else {
	i = _hsize ;
      }
    }

    return ;
  }

}

/** @} */ //end of group class.

#endif	/* PRIORITY_QUEUE_HPP */
