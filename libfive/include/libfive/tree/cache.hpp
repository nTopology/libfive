/*
libfive: a CAD kernel for modeling with implicit functions

Copyright (C) 2017  Matt Keeter

This Source Code Form is subject to the terms of the Mozilla Public
License, v. 2.0. If a copy of the MPL was not distributed with this file,
You can obtain one at http://mozilla.org/MPL/2.0/.
*/
#pragma once

#include <map>
#include <memory>
#include <mutex>
#include <assert.h>

#include "libfive/export.hpp"
#include "libfive/tree/tree.hpp"

namespace libfive {

/*
 *  A Cache stores values in a deduplicated math expression.  It is 
 *  static-initialization-safe: If its methods are called before it is 
 *  initialized or after deinitialization (e.g. due to destroying a static
 *  std::shared_ptr<Tree> that was initialized before the cache), the this 
 *  pointer is not used, and instead a non-deduplicated tree is returned if 
 *  necessary; this also occurs if accessed from a Handle that could not 
 *  lock the mutex due to concerns that it may not have been initialized, 
 *  or if the methods are called from Cache*(nullptr).
 */
class Cache
{
    /*  Helper typedef to avoid writing this over and over again  */
    typedef std::shared_ptr<Tree::Tree_> Node;

    /*  Handle to safely access cache  */
    class Handle
    {
    public:
        Handle() 
        {
          if (_exists) {
            lock = std::unique_lock<std::recursive_mutex>(mut);
          }
          // Otherwise, mut might not exist yet or may have already
          // been destroyed, so accessing it may not be safe.
        }
        Cache* operator->() const 
        { 
          return lock.owns_lock() ? &_instance : nullptr; 
        }
    protected:
        std::unique_lock<std::recursive_mutex> lock;
    };

public:
    /*
     *  Returns a safe (locking) handle to the global Cache
     */
    static Handle instance() { return Handle(); }

    Node constant(float v);
    Node operation(Opcode::Opcode op, Node lhs=nullptr, Node rhs=nullptr,
                   bool simplify=true);

    Node X() { return operation(Opcode::VAR_X); }
    Node Y() { return operation(Opcode::VAR_Y); }
    Node Z() { return operation(Opcode::VAR_Z); }

    Node var();

    /*
     *  Called when the last Tree_ is destroyed
     */
    void del(float v);
    void del(Opcode::Opcode op, Node lhs=nullptr, Node rhs=nullptr);

    /*
     *  Returns the given node as an affine sum-of-multiplications
     *
     *  This is a building block for automatic collapsing of affine
     *  expressions, exposed here primarily for unit testing.
     */
    std::map<Node, float> asAffine(Node n);

    /*
     *  Converts a sum-of-multiplications into an affine tree.
     */
    Node fromAffine(const std::map<Node, float>& ns);

protected:
    /*
     *  Cache constructor is private so outsiders must use instance()
     */
    Cache() 
    {
      std::unique_lock<std::recursive_mutex> lock(mut);
      _exists = true; 
    }

    ~Cache() 
    { 
      std::unique_lock<std::recursive_mutex> lock(mut);
      _exists = false; 
    }

    /*
     *  Checks if this is a valid cache.  May be safely called on any
     *  pointer.
     */
    bool isValid() {
      assert(this == nullptr || this == &_instance);
      return this != nullptr && _exists;
    }

    /*
     *  Checks whether the operation is an identity operation
     *  If so returns an appropriately simplified tree
     *  i.e. (X + 0) will return X
     */
    Node checkIdentity(Opcode::Opcode op, Node a, Node b);

    /*
     *  If the opcode is commutative, consider tweaking tree structure
     *  to keep it as balanced as possible.
     */
    Node checkCommutative(Opcode::Opcode op, Node a, Node b);

    /*
     *  Checks whether a `op` b can be replaced by a simpler affine
     *  form, e.g. (2*x + y) - (4*y) --> 2*x - 3*y
     */
    Node checkAffine(Opcode::Opcode op, Node a, Node b);

    /*
     *  A Key uniquely identifies an operation Node, so that we can
     *  deduplicate based on opcode  and arguments
     */
    typedef std::tuple<Opcode::Opcode,  /* opcode */
                       Tree::Id,        /* lhs */
                       Tree::Id         /* rhs */ > Key;
    std::map<Key, std::weak_ptr<Tree::Tree_>> ops;

    /*  Constants in the tree are uniquely identified by their value  */
    std::map<float, std::weak_ptr<Tree::Tree_>> constants;

    /*  nan cannot be stored in the usual map, so the nan constant lives here */
    std::weak_ptr<Tree::Tree_> nan_constant;

    /*  Oracles do not need to use the cache to be deduplicated, since they
     *  are created from unique_ptr's, and therefore are already impossible
     *  to duplicate.  */

    static FIVE_EXPORT std::recursive_mutex mut;
    static FIVE_EXPORT bool _exists;
    static FIVE_EXPORT Cache _instance;
};

}   // namespace libfive
