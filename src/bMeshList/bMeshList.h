/***************************************************************************
**
**  		     tRIBS Distributed Hydrologic Model
**
**              TIN-based Real-time Integrated Basin Simulator
**		         Ralph M. Parsons Laboratory
**  		    Massachusetts Institute of Technology
**  
**
**  bMeshList.h: Header file for derived classes bMeshList and bMeshListIter
**
**  A bMeshList is derived from the generic linked list class tList.
**  It is used in CHILD to store lists of grid elements (nodes and edges),
**  and differs from a generic list in being divided into two parts:
**  (1) an "active" part, representing elements that are not part of the
**  mesh boundary and are therefore subject to active processes (whatever
**  those may be; in CHILD the processes are runoff, erosion, and
**  sedimentation); and (2) a "boundary" part, containing elements along
**  the mesh boundary.
**
**  A bMeshListIter is an iterator for a bMeshList. It has the same services
**  as a tListIter. It also will move to the last "active" (non-boundary) 
**  node on a grid list, or to the first boundary node on the list. It adds 
**  special functions FirstP, NextP that are identical to the tListIter 
**  functions First and Next except that they return a pointer to the data
**  portion of the node (or zero if the end of the list is reached, or the 
**  current node is null for some other reason).
**
**************************************************************************/

#ifndef BMESHLIST_H
#define BMESHLIST_H

#include "src/Headers/Inclusions.h"
#include <list>

//=========================================================================
//
//
//                  Section 1: bMeshList Class Declarations
//
//
//=========================================================================

/**************************************************************************
**
** bMeshList()
**
** Class bMeshList implements a linked list that is divided into two
** parts, an "active" (front) and "inactive" (back) part. Derived from tList.
**
**************************************************************************/

template<class NodeType>
class bMeshList : public std::list<NodeType>
{
public:
  bMeshList();
  ~bMeshList();

  int getActiveSize() const	{ return nActiveNodes; }
  int getTotalSize() const	{ return totalNodes; }

  std::_List_iterator<NodeType> getLastActive() const;

  int isActiveEmpty() const;
  int isBoundEmpty() const;

  void insertAtBack(const NodeType);
  void insertAtBoundFront(const NodeType);
  int removeFromBoundFront(NodeType&);

  void insertAtActiveBack(const NodeType);
  int removeFromActiveBack(NodeType&);

  int removeActive(std::_List_iterator<NodeType>);
  int removeNext(NodeType &value, std::_List_iterator<NodeType>);
  int removePrev(NodeType &value, std::_List_iterator<NodeType>);

  void moveToFront(std::_List_iterator<NodeType>);
  void moveToBack(std::_List_iterator<NodeType>);
  void moveToBack(NodeType);

  void moveToActiveBack(std::_List_iterator<NodeType>);
  void moveToBoundFront(std::_List_iterator<NodeType>);

  void insertAtFront(const NodeType);
  int removeFromFront(NodeType&);

  int InActiveList(std::_List_iterator<NodeType>);
  int InActiveList(NodeType);

  void Flush();
   
protected:
  int nActiveNodes;
  int totalNodes;
  std::_List_iterator<NodeType> lastactive;
};

//=========================================================================
//
//
//                  Section 1: bMeshList Constructors/Destructors
//
//
//=========================================================================
// defined in bMeshList.h
template< class NodeType >
bMeshList<NodeType>::bMeshList() :
        nActiveNodes(0),
        totalNodes(0),
        lastactive(this->end())
{}

template< class NodeType >
bMeshList<NodeType>::~bMeshList(){}

//=========================================================================
//
//
//                  Section 2: bMeshList Functions
//
//
//=========================================================================

template< class NodeType >
std::_List_iterator<NodeType> bMeshList<NodeType>::getLastActive() const
{
    std::_List_iterator<NodeType> firstboundary = lastactive;
    firstboundary++;
    return firstboundary;
}

template< class NodeType >
int bMeshList<NodeType>::isActiveEmpty() const
{
    if (nActiveNodes == 0)
        return 1;
    else
        return 0;
}

template< class NodeType >
int bMeshList<NodeType>::isBoundEmpty() const
{
    if (totalNodes - nActiveNodes == 0)
        return 1;
    else
        return 0;
}


/**************************************************************************
**
**  bMeshList insertion and removal functions
**
**  Adds and removes items to/from the list. Supplements std::list
**  functionality by adding capability to add items to front of
**  "boundary" section or rear of "active" section. Updates
**  nActiveNodes as appropriate.
**
**************************************************************************/

template< class NodeType >
void bMeshList<NodeType>::insertAtBack(const NodeType value)
{
    this->push_back(value);
    totalNodes++;
}

template< class NodeType >
void bMeshList<NodeType>::insertAtFront(const NodeType value)
{
    this->insert(this->begin(), value);
    if (nActiveNodes == 0) {
        lastactive = this->begin();
    }
    nActiveNodes++;
    totalNodes++;
}

template< class NodeType >
int bMeshList<NodeType>::removeFromFront(NodeType& value)
{
    // Active is empty and boundary is empty
    if (this->empty())
        return 0;

    // If the node at the front is active and not boundary
    if (nActiveNodes != 0) {
        nActiveNodes--;
        totalNodes--;
        if (nActiveNodes == 0)
            lastactive = this->end();
    }

    // Remove the node at front either active or boundary
    std::_List_iterator<NodeType> nodeToRemove = this->begin();
    value = (*nodeToRemove);
    this->erase(nodeToRemove);
    totalNodes--;
    return 1;
}


template< class NodeType >
void bMeshList<NodeType>::insertAtBoundFront(const NodeType value)
{
    assert( this != 0 );

    // Active is empty and boundary is empty
    if (this->empty()) {
        this->push_back(value);
        totalNodes++;
    }

        // Active is empty and boundary is not empty
    else if (lastactive == this->end()) {
        this->insert(this->begin(), value);
        totalNodes++;
    }

        // Active is not empty
    else {
        std::_List_iterator<NodeType> firstbound = lastactive;
        firstbound++;
        this->insert(firstbound, value);
        totalNodes++;
    }
}

template< class NodeType >
int bMeshList<NodeType>::removeFromBoundFront(NodeType &value)
{
    assert(&value != 0);

    // Active and boundary are both empty
    if (this->empty())
        return 0;

        // Active is not empty and boundary is empty
    else if (totalNodes == nActiveNodes)
        return 0;

    // Active is empty or not and boundary is not empty
    std::_List_iterator<NodeType> firstbound;
    if (nActiveNodes == 0)
        firstbound = this->begin();
    else {
        firstbound = lastactive;
        firstbound++;
    }

    value = (*firstbound);
    this->erase(firstbound);
    totalNodes--;
    return 1;
}

template< class NodeType >
void bMeshList<NodeType>::insertAtActiveBack(const NodeType value)
{
    assert( this != 0 );

    // Active and boundary are both empty
    if (this->empty()) {
        this->push_back(value);
        lastactive = this->begin();
    }

        // Active is empty and boundary is not empty
    else if (nActiveNodes == 0) {
        this->insert(this->begin(), value);
        lastactive = this->begin();
    }

        // Active is not empty and boundary is empty
    else if (totalNodes == nActiveNodes) {
        this->push_back(value);
        lastactive++;
    }

        // Active is not empty and boundary is not empty
    else {
        std::_List_iterator<NodeType> firstboundary = lastactive;
        firstboundary++;
        this->insert(firstboundary, value);
        lastactive++;
    }
    nActiveNodes++;
    totalNodes++;
}

template< class NodeType >
int bMeshList<NodeType>::removeFromActiveBack(NodeType &value)
{
    // Active is empty
    if (nActiveNodes == 0)
        return 0;

    // Active has one node
    if (nActiveNodes == 1) {
        value = (*lastactive);
        this->erase(lastactive);
        lastactive = this->end();
    }

        // Active has more than one node
    else {
        value = (*lastactive);
        std::_List_iterator<NodeType> nodeToRemove = lastactive;
        lastactive--;
        this->erase(nodeToRemove);
    }
    nActiveNodes--;
    totalNodes--;
    return 1;
}

// delete active node under the iterator
template< class NodeType >
int bMeshList<NodeType>::removeActive(std::_List_iterator<NodeType> ptr)
{
    // If the passed iterator is empty
    if (ptr == this->end())
        return 0;

    // If node to remove is the last active
    if (ptr == lastactive)
        return removeFromActiveBack((*ptr));

    // Node to remove is in middle of active
    this->erase(ptr);
    nActiveNodes--;
    totalNodes--;
    return 1;
}

//delete next node
template< class NodeType >
int bMeshList<NodeType>::removeNext(NodeType &value,
                                    std::_List_iterator<NodeType> ptr)
{
    // If the passed iterator is empty
    if (ptr == this->end())
        return 0;

    // If the iterator after the passed iterator is empty
    std::_List_iterator<NodeType> nodeToRemove = ptr;
    nodeToRemove++;
    if (nodeToRemove == this->end())
        return 0;

    // If node to remove is the last active
    if (nodeToRemove == lastactive)
        return removeFromActiveBack(value);

    // If node to remove is the first boundary
    std::_List_iterator<NodeType> temp = lastactive;
    temp++;
    if (nodeToRemove == temp)
        return removeFromBoundFront(value);

    // Node to remove is in middle of active or middle of boundary
    value = (*nodeToRemove);
    if (value.getBoundaryFlag() == 0)
        nActiveNodes--;
    this->erase(nodeToRemove);
    totalNodes--;
    return 1;
}

//delete previous node
template< class NodeType >
int bMeshList<NodeType>::removePrev(NodeType &value,
                                    std::_List_iterator<NodeType> ptr)
{
    // If the passed iterator is empty
    if (ptr == this->end())
        return 0;

    // If the passed iterator is first
    if (ptr == this->begin())
        return 0;

    // If the passed iterator is the first boundary
    std::_List_iterator<NodeType> nodeToRemove = ptr - 1;
    if (nodeToRemove == lastactive)
        return removeFromActiveBack(value);

    // Node to remove is in middle of active or middle of boundary
    value = (*nodeToRemove);
    if (value.getBoundaryFlag() == 0)
        nActiveNodes--;
    this->erase(nodeToRemove);
    totalNodes--;
    return 1;
}

/**************************************************************************
**
**  bMeshList::moveToBack ( tListNode * )
**
**  Moves niter to the back of the list (the boundary portion).
**  Handles case of moved node being the last active node, in which case
**  _lastactive_ needs to be updated.
**
**************************************************************************/

template< class NodeType >
void bMeshList<NodeType>::moveToBack(std::_List_iterator<NodeType> niter)
{
    assert(niter != this->end());

    // Node to move is active
    if (InActiveList(niter)) {
        nActiveNodes--;
        // Node to move is the last active
        if (niter == lastactive) {
            // Node to move is the only active
            if (nActiveNodes == 0)
                lastactive = this->end();
            else
                lastactive--;
        }
    }

    // Move the node to the end of boundary
    NodeType value = (*niter);
    this->erase(niter);
    this->push_back(value);
}


/**************************************************************************
**
**  bMeshList::moveToBack ( NodeType )
**
**  Finds the ListNode whose data are identical to node and calls
**  moveToBack( tListNode ) to move it to the back of the list.
**
**************************************************************************/

template< class NodeType >
void bMeshList<NodeType>::moveToBack(NodeType node)
{
    std::_List_iterator<NodeType> nodeToMove =
            find(this->begin(), this->end(), node);

    if (nodeToMove != this->end())
        moveToBack(nodeToMove);
}

/**************************************************************************
**
**  bMeshList::moveToFront()
**
**  Moves niter to the front of the list, taking care to handle the case
**  in which the node being moved is the last on the active section
**
**************************************************************************/

template< class NodeType >
void bMeshList<NodeType>::moveToFront(std::_List_iterator<NodeType> niter)
{
    // Node to move is not already in the front
    if (niter == this->begin())
        return;

    // Node to move is the last active
    if (niter == lastactive)
        lastactive--;

    // Node to move is boundary
    if (!InActiveList(niter))
        nActiveNodes++;

    NodeType value = (*niter);
    this->erase(niter);
    this->insert(this->begin(), value);
}


/**************************************************************************
**
**  bMeshList::moveToActiveBack()
**
**  Moves niter to the back of the "active" portion of the list
**  (does not update nActiveNodes if the node happens to be inactive!)
**
**************************************************************************/

template< class NodeType >
void bMeshList<NodeType>::moveToActiveBack(std::_List_iterator<NodeType> niter)
{
    // Node to move is already at active back
    if (niter == lastactive)
        return;

    // Node to move is active
    if (InActiveList(niter)) {
        NodeType value = (*niter);
        insertAtActiveBack(value);	// will increment nActiveNodes
        this->erase(niter);
        nActiveNodes--;
    }

        // Node to move is boundary
    else {
        NodeType value = (*niter);
        insertAtActiveBack(value);
        this->erase(niter);
    }
}

/**************************************************************************
**
**  bMeshList::moveToBoundFront()
**
**  Moves niter to the front of the "boundary" portion of the list,
**  making sure to update nActiveNodes is the node was previously on
**  the active portion of the list.
**
**************************************************************************/

template< class NodeType >
void bMeshList<NodeType>::moveToBoundFront(std::_List_iterator<NodeType> niter)
{
    // Node to move is already at bound front
    std::_List_iterator<NodeType> temp = lastactive;
    temp++;
    if (niter == temp)
        return;

    // Node to move is last active
    if (niter == lastactive) {
        lastactive--;
        nActiveNodes--;
        return;
    }

    // Node to move is active
    NodeType value = (*niter);
    if (InActiveList(niter)) {
        nActiveNodes--;
        this->insert(lastactive, value);
        this->erase(niter);
    }

        // Node to move is boundary
    else {
        insertAtBoundFront(value);
        this->erase(niter);
    }
}


/**************************************************************************
**
**  bMeshList::Flush()
**
**  Also reinitializes lastactive and nActiveNodes
**
**************************************************************************/

template< class NodeType >
void bMeshList<NodeType>::Flush()
{
    this->erase(this->begin(), this->end());
    lastactive = this->end();
    nActiveNodes = 0;
}

/**************************************************************************
**
**  bMeshList::InActiveList
**
**  Reports whether a given list node is in the active portion of the list.
**
**  Parameters:  niter -- list node to test
**  Returns:  1 if niter is present in the active portion of the list,
**            0 otherwise.
**
**************************************************************************/

template< class NodeType >
int bMeshList<NodeType>::InActiveList(std::_List_iterator<NodeType> niter )
{
    if (nActiveNodes == 0)
        return 0;

    std::_List_iterator<NodeType> iter = this->begin();
    while (iter != this->getLastActive() && iter != niter)
        iter++;

    if (iter == niter)
        return 1;
    else
        return 0;
}

template< class NodeType >
int bMeshList<NodeType>::InActiveList(NodeType niter)
{
    if (nActiveNodes == 0)
        return 0;

    std::_List_iterator<NodeType> iter = this->begin();
    while (iter != this->getLastActive() && (*iter) != niter)
        iter++;

    if ((*iter) == niter)
        return 1;
    else
        return 0;
}


#endif

//=========================================================================
//
//
//                      End of bMeshList.h
//
//
//=========================================================================
