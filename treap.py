# This module implements a simple balanced binary search tree (Cartesian tree aka Treap)
# It supports insertion and query for the number and sum of elements which have value
# lower than the average value of inserted items.

import random

class TreapNode:
    '''
    Class representing a node of a tree
    '''
    def __init__(self):
        self.element = 0.0      # value of stored element
        self.subtree_sum = 0.0  # sum of values in a subtree (including this node)
        self.subtree_cnt = 0    # number of nodes in a subtree (including this node)
        self.weight = 0.0       # balancing weight
        self.left_node = -1     # index of left child (-1 if there are none)
        self.right_node = -1    # index of right child (-1 if there are none)

class Treap:
    '''
    Class representing a balanced binary search tree
    '''
    def __init__(self):
        self.nodes = []     # storage of tree nodes
        self.root = -1      # index of root node (-1 means no root)
    
    def next_weight():
        '''
        Generates a new uniformly distributed random weight
        '''
        return random.randint(-1e6, 1e6)
    
    def update_subtree(self, node):
        '''
        Updates subtree_sum and subtree_cnt of a node with index 'node' assuming that
        all children have valid values, O(1)
        '''
        self.nodes[node].subtree_sum = self.nodes[node].element
        self.nodes[node].subtree_cnt = 1
        if self.nodes[node].left_node != -1:
            self.nodes[node].subtree_sum += self.nodes[self.nodes[node].left_node].subtree_sum
            self.nodes[node].subtree_cnt += self.nodes[self.nodes[node].left_node].subtree_cnt
        if self.nodes[node].right_node != -1:
            self.nodes[node].subtree_sum += self.nodes[self.nodes[node].right_node].subtree_sum    
            self.nodes[node].subtree_cnt += self.nodes[self.nodes[node].right_node].subtree_cnt    
    
    def split(self, node, element):
        '''
        Splits the subtree with root at 'node' into two unconnected trees, the first containing
        all elements lower than 'element'. Returns two roots, O(log n)
        '''
        if node == -1:
            return -1, -1
        else:
            if element < self.nodes[node].element:
                l, r = self.split(self.nodes[node].left_node, element)
                self.nodes[node].left_node = r
                self.update_subtree(node)
                return l, node
            else:
                l, r = self.split(self.nodes[node].right_node, element)
                self.nodes[node].right_node = l
                self.update_subtree(node)
                return node, r
    
    def join(self, left, right):
        '''
        Joins two trees with root indices at 'left' and 'right' and returns a single root, O(log n)
        '''
        if left == -1 or right == -1:
            return left if left != -1 else right
        else:
            if self.nodes[left].weight > self.nodes[right].weight:
                self.nodes[left].right_node = self.join(self.nodes[left].right_node, right)
                self.update_subtree(left)
                return left
            else:
                self.nodes[right].left_node = self.join(left, self.nodes[right].left_node)
                self.update_subtree(right)
                return right
    
    def add_element(self, element):
        '''
        Adds a value 'element' in the tree (creating a new node), O(log n)
        '''
        new_node = TreapNode()
        new_node.element = element
        new_node.weight = Treap.next_weight()
        self.nodes.append(new_node)
        self.root = self.add_element_internal(self.root)
    
    def add_element_internal(self, node):
        '''
        Adds a node inside a subtree with root 'node', assuming that inserted element 
        is positioned at the back of self.nodes list, O(log n)
        '''
        if node == -1:
            self.update_subtree(-1)
            return len(self.nodes) - 1
        else:
            if self.nodes[-1].weight > self.nodes[node].weight:
                l, r = self.split(node, self.nodes[-1].element)
                self.nodes[-1].left_node = l
                self.nodes[-1].right_node = r
                self.update_subtree(-1)
                return len(self.nodes) - 1
            else:
                if self.nodes[-1].element < self.nodes[node].element:
                    self.nodes[node].left_node = self.add_element_internal(self.nodes[node].left_node)
                else:
                    self.nodes[node].right_node = self.add_element_internal(self.nodes[node].right_node)
                self.update_subtree(node)
                return node
    
    def get_lower(self, node, value):
        '''
        Gets number and sum of elements in the tree with root 'node' which have value less
        than 'value', O(log n)
        '''
        if node == -1:
            return 0, 0
        if self.nodes[node].element < value:
            c, s = self.get_lower(self.nodes[node].right_node, value)
            if self.nodes[node].left_node != -1:
                c += self.nodes[self.nodes[node].left_node].subtree_cnt
                s += self.nodes[self.nodes[node].left_node].subtree_sum
            return c + 1, s + self.nodes[node].element
        else:
            return self.get_lower(self.nodes[node].left_node, value)