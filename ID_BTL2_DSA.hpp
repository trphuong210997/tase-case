#include "main.hpp"
// #include "Dataset.hpp"
/* TODO: Please design your data structure carefully so that you can work with the given dataset
 *       in this assignment. The below structures are just some suggestions.
 */
struct Node
{
    vector<int> data;
    Node *left;
    Node *right;
    Node(vector<int> data, Node *left = nullptr, Node *right = nullptr)
    {
        this->data = data;
        this->left = left;
        this->right = right;
    }

    void print() const
    {
        OUTPUT << "(";
        for(int i = 0; i < data.size(); i++)
        {
            OUTPUT << data[i];
            if (i == data.size() - 1) 
                OUTPUT << ")";
            else 
                OUTPUT << ",";
        }
    }
};

class kDTree
{
private:
    int count;
    int k;
    Node *root;

public:
    kDTree(int k = 2) : k(k), count(0), root(nullptr) {}

    kDTree(int k, int count, Node *root) : k(k), count(count), root(root) {}

    ~kDTree() {
        delete root;
    }

    Node *copyData(Node *node) {
        if (node == nullptr) 
            return nullptr;
        
        return new Node(node->data, copyData(node->left), copyData(node->right));
    }

    const kDTree &operator=(const kDTree &other) {
        this->k = other.k;
        this->count = other.count;
        this->root = copyData(other.root);
        return *this;
    }

    kDTree(const kDTree &other) {
        this->k = other.k;
        this->count = other.count;
        this->root = copyData(other.root);
    }

    void LNR(Node *node) const {
        if (node == nullptr)
            return;
        
        LNR(node->left);
        node->print();
        OUTPUT << " ";
        LNR(node->right);
    }

    void NLR(Node *node) const {
        if (node == nullptr)
            return;
        
        node->print();
        // if (node->left != nullptr)
            OUTPUT << " ";
        NLR(node->left);
        // if (node->right != nullptr)
        //     OUTPUT << " ";
        NLR(node->right);
    }

    void LRN(Node *node) const {
        if (node == nullptr)
            return;
        
        LRN(node->left);
        // if (node->left != nullptr)
        //     OUTPUT << " ";
        LRN(node->right);
        // if (node->right != nullptr)
            // OUTPUT << " ";
        node->print();
        OUTPUT << " ";
    }

    void inorderTraversal() const {
        return LNR(root);
    }

    void preorderTraversal() const {
        return NLR(root);
    }

    void postorderTraversal() const {
        return LRN(root);
    }

    int calHeight(Node *node) const {
        if (node == nullptr) 
            return 0;

        if (node->left == nullptr && node->right == nullptr)
            return 1;
        
        return 1 + max(calHeight(node->left), calHeight(node->right));
    }

    int height() const {
        return calHeight(root);
    }

    int nodeCount() const {
        return count;
    }

    int TreeLeafCount(Node *node) const {
        if (node == nullptr)
            return 0;

        if (node != nullptr && node->left == nullptr && node->right == nullptr)
            return 1;
        
        return TreeLeafCount(node->left) + TreeLeafCount(node->right);
    }

    int leafCount() const {
        return TreeLeafCount(root);
    }

    Node *insertHelper(Node *node, const vector<int> &point, int treeLevel) {
        if (node == nullptr)
            return new Node(point);
        
        if (point[treeLevel % k] < node->data[treeLevel % k])
            node->left = insertHelper(node->left, point, treeLevel + 1);
        else
            node->right = insertHelper(node->right, point, treeLevel + 1);
        
        return node;
    }
    

    void insert(const vector<int> &point) {
        this->root = insertHelper(root, point, 0);
        this->count++;
    }

    bool searchHelper(Node *node, const vector<int> &point, int treeLevel) {
        if (node == nullptr)
            return false;
        
        if (node->data == point)
            return true;

        if (point[treeLevel % k] < node->data[treeLevel % k])
            return searchHelper(node->left, point, treeLevel + 1);
        else
            return searchHelper(node->right, point, treeLevel + 1);
    }

    bool search(const vector<int> &point) {
        return searchHelper(root, point, 0);
    }

    Node *minNode(Node *x, Node *y, Node *z, int alpha) {
        Node *min = x;
        if (y != nullptr && y->data[alpha] < min->data[alpha])
            min = y;
        if (z != nullptr && z->data[alpha] < min->data[alpha])
            min = z;
        return min;
    }

    Node *findMin(Node *node, int treeLevel, int alpha) {
        if (node == nullptr)
            return node;
        
        int dimension = treeLevel % k;

        if (dimension == alpha) {
            if (node->left != nullptr) {
                return findMin(node->left, treeLevel + 1, alpha);
            }
            else
                return node;
        }

        else 
            return minNode(node, findMin(node->left, treeLevel + 1, alpha), findMin(node->right, treeLevel + 1, alpha), alpha);
    }

    void removeHelper(Node* &node, Node* prev_node, const vector<int> &point, int treeLevel) {
        if (node == nullptr)
            return;
        if (node->data == point) {
            int dimension = treeLevel % k;

            if (node->left == nullptr && node->right == nullptr) {
                if (prev_node != nullptr) {
                    if (prev_node->left == node)
                        prev_node->left = nullptr;
                    else if (prev_node->right == node)
                        prev_node->right = nullptr;
                }
                delete node;
                node = nullptr;
                prev_node = nullptr;
                count--;
                return;
            }

            if (node->right != nullptr) {
                Node *minNodeRight = findMin(node->right, treeLevel + 1, dimension);
                node->data = minNodeRight->data;
                return removeHelper(node->right, node, minNodeRight->data, treeLevel + 1); 
            }

            if (node->right == nullptr) {
                Node *minNodeLeft = findMin(node->left, treeLevel + 1, dimension);
                node->data = minNodeLeft->data;
                node->right = node->left;
                node->left = nullptr;
                return removeHelper(node->right, node, minNodeLeft->data, treeLevel + 1);
            }
        }

        if (point[treeLevel % k] < node->data[treeLevel % k])
            return removeHelper(node->left, node, point, treeLevel + 1);
        else
            return removeHelper(node->right, node, point, treeLevel + 1);
    }

    void remove(const vector<int> &point) {
        return removeHelper(root, nullptr, point, 0);
    }

    void sort(vector<vector<int>> &pointList, int k, int dimension) {
        for (int i = 0; i < pointList.size(); i++) {
            for (int j = i + 1; j < pointList.size(); j++) {
                if (pointList[i][dimension] > pointList[j][dimension]) {
                    vector<int> temp = pointList[i];
                    pointList[i] = pointList[j];
                    pointList[j] = temp;
                }
            }
        }
    }

    int findFirstMedian (const vector<vector<int>> &pointList) {
        int ListMidpoint;
        if (pointList.size() % 2 == 0)
            ListMidpoint = pointList.size() / 2 - 1;
        else
            ListMidpoint = pointList.size() / 2;
        for (int i = 0; i < ListMidpoint; i++) {
            if (pointList[i] == pointList[ListMidpoint])
                return i;
        }
        return ListMidpoint;
    }
    
    Node *buildTreeHelper(const vector<vector<int>> &pointList, int treeLevel) {
        vector<vector<int>> pointListCopy = pointList;
        if (pointList.size() == 0)
            return nullptr;
        
        int dimension = treeLevel % k;
        sort(pointListCopy, k, dimension);
        int ListMidpoint = findFirstMedian(pointListCopy);
        Node *node = new Node(pointListCopy[ListMidpoint]);
        count++;

        node->left = buildTreeHelper(vector<vector<int>>(pointListCopy.begin(), pointListCopy.begin() + ListMidpoint), treeLevel + 1);
        node->right = buildTreeHelper(vector<vector<int>>(pointListCopy.begin() + ListMidpoint + 1, pointListCopy.end()), treeLevel + 1);

        return node;
    }

    void buildTree(const vector<vector<int>> &pointList) {
        this->root = buildTreeHelper(pointList, 0);
    }

    double calculateDistance(const vector<int> &point1, const vector<int> &point2) {
        double distance = 0;
        for (int i = 0; i < point1.size(); i++) {
            distance += pow(point1[i] - point2[i], 2);
        }
        return sqrt(distance);
    }

    // Node *findBestNode(Node *node, const vector<int> &target, int treeLevel) {
    //     if (node)
    //         return nullptr;
        
    //     if (target[treeLevel % k] < node->data[treeLevel % k]) 
    //         return findBestNode(node->left, target, treeLevel + 1);
    //     else
    //         return findBestNode(node->right, target, treeLevel + 1);
        
    //     return best;
    // }
    
    // void nearestNeighbourHelper(Node *node, const vector<int> &target, Node *best, int treeLevel) {
    //     Node *bestNode = findBestNode(root, target, 0);
    // }

    // void nearestNeighbour(const vector<int> &target, Node *best) {
    //     best = nullptr;
    //     return nearestNeighbourHelper(root, target, best, 0);
    // }

    void kNearestNeighbour(const vector<int> &target, int k, vector<Node *> &bestList);
};



// Please add more or modify as needed
