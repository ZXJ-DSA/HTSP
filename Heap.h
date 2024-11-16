// Basic classes, Timer, Semaphore, Heap
// created by Mengxuan, modified by Xinjie

#pragma once

#include <mutex>
#include <condition_variable>
#include <vector>
#include <chrono>

struct Timer
{
private:
    std::chrono::high_resolution_clock::time_point t1, t2;//varibles for time record
//    std::chrono::steady_clock::time_point t1, t2;//varibles for time record
    std::chrono::duration<double> time_span;
//    std::chrono::time_point<std::chrono::high_resolution_clock> m_begin, m_end;
public:
    void start()
    {
        t1 = std::chrono::high_resolution_clock::now();
//        t1 = std::chrono::steady_clock::now();
//        m_begin = std::chrono::high_resolution_clock::now();
    }
    void stop()
    {
        t2 = std::chrono::high_resolution_clock::now();
//        t2 = std::chrono::steady_clock::now();
//        m_end = std::chrono::high_resolution_clock::now();
    }
    double GetRuntime()//return time in second
    {
        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);//std::chrono::nanoseconds
        return time_span.count();
    }
//    int64_t GetRuntime2() const//return time in second
//    {
//        time_span = std::chrono::duration_cast<std::chrono::duration<double>>(t2 - t1);//std::chrono::nanoseconds
//        return time_span.count();
//        return std::chrono::duration_cast<std::chrono::nanoseconds>(m_end - m_begin).count();//
//    }
};

class Semaphore
{
public:
    Semaphore (int count_ = 0)
            : count(count_) {}

    inline void notify()
    {
        std::unique_lock<std::mutex> lock(mtx);
        count++;
        cv.notify_one();
    }

    inline void wait()
    {
        std::unique_lock<std::mutex> lock(mtx);

        while(count == 0)
        {
            cv.wait(lock);

        }
        count--;
    }

private:
    std::mutex mtx;
    std::condition_variable cv;
    int count;//the number of thread waiting here?
};

namespace benchmark {

#define NULLINDEX 0xFFFFFFFF

    template<int log_k, typename id_t, typename k_t >//id,value
    class heap {//binary heap
    public:

        // Expose types.
        typedef k_t key_t;
        typedef id_t node_t;

        // Some constants regarding the elements.
        //static const node_t NULLINDEX = 0xFFFFFFFF;
        static const node_t k = 1 << log_k;//equals k = 1*2^log_k, usually log_k = 2, so k = 4

        // A struct defining a heap element.
        struct element_t {
            key_t key;
            node_t element;

            element_t() : key(0), element(0) {}

            element_t(const key_t k, const node_t e) : key(k), element(e) {}
        };


        //public:

        // Constructor of the heap.
        heap(node_t n) : n(0), max_n(n), elements(n), position(n, NULLINDEX) {//n is the number of elements in current
            // state, max_n is the size of heap
        }

        heap(): n(0), max_n(0), elements(0), position(0, NULLINDEX) {}

        ~heap(){}

        // Risize the heap
        inline void resize(node_t a){
            n = 0; max_n = a;
            elements.resize(a);
            position.resize(a, NULLINDEX);
        }

        // Size of the heap.
        inline node_t size() const {
            return n;
        }

        // Heap empty?
        inline bool empty() const {
            return size() == 0;
        }

        // Extract min element.
        inline void extract_min(node_t &element, key_t &key) {
            assert(!empty());

            element_t &front = elements[0];

            // Assign element and key.
            element = front.element;
            key = front.key;

            // Replace elements[0] by last element.
            position[element] = NULLINDEX;
            --n;//n=n-1
            if (!empty()) {
                front = elements[n];//elements[n] is the top element
                position[front.element] = 0;//make its position valid, it is also the smallest one
                sift_down(0);
            }
        }

        inline key_t top_key() {//get the key, i.e. minimal cost
            assert(!empty());

            element_t &front = elements[0];

            return front.key;

        }

        inline node_t top_id() {//get the value, i.e. id number of minimal cost

            assert(!empty());

            element_t &front = elements[0];

            return front.element;
        }

        // Update an element of the heap.
        inline void update(const node_t element, const key_t key) {

            if (position[element] == NULLINDEX) {//if originally NULL
                element_t &back = elements[n];//add new element to position n
                back.key = key;
                back.element = element;
                position[element] = n;//set position id to n
                sift_up(n++);
            } else {//if already valid, update the value
                node_t el_pos = position[element];//position information
                element_t &el = elements[el_pos];//get the element
                if (key > el.key) {//update the element
//                if (key > el.key || (key <= el.key && element > el.element)) {//update the element || (elements[parent_i].key <= elements[cur_i].key && elements[parent_i].element > elements[cur_i].element)
                    el.key = key;
                    sift_down(el_pos);
                } else {
                    el.key = key;
                    sift_up(el_pos);
                }
            }
        }

        // Clear the heap.
        inline void clear() {
            for (node_t i = 0; i < n; ++i) {
                position[elements[i].element] = NULLINDEX;
            }
            n = 0;
        }

        // Cheaper erase.
        inline void erase(node_t v) {
            position[v] = NULLINDEX;
        }

        inline void clear_n() {
            n = 0;
        }

        // Test whether an element is contained in the heap.
        inline bool contains(const node_t element) const {
            return position[element] != NULLINDEX;
        }

        //return current elements information
        void get_elements(std::vector<std::pair<int,int>> &e_vector){
            std::pair<int,int> temp_pair;

            for(int i=0;i<n;i++){
                temp_pair.first = elements[i].key;
                temp_pair.second = elements[i].element;
                e_vector.push_back(temp_pair);
            }
        }

    protected:

        // Sift up an element.
        inline void sift_up(node_t i) {
            assert(i < n);
            node_t cur_i = i;
            while (cur_i > 0) {
                node_t parent_i = (cur_i - 1) >> log_k;//equals (cur_i - 1)/(2^log_k)
                if (elements[parent_i].key > elements[cur_i].key)//compare with parent node, if smaller, then swap
//                if (elements[parent_i].key > elements[cur_i].key || (elements[parent_i].key <= elements[cur_i].key && elements[parent_i].element > elements[cur_i].element))//compare with parent node, if smaller, then swap
                    swap(cur_i, parent_i);
                else
                    break;
                cur_i = parent_i;
            }
        }

        // Sift down an element.
        inline void sift_down(node_t i) {
            assert(i < n);

            while (true) {
                node_t min_ind = i;
                key_t min_key = elements[i].key;

                node_t child_ind_l = (i << log_k) + 1;//equals i*2^log_k + 1
                node_t child_ind_u = std::min(child_ind_l + k, n);//equals min(child_ind_l+4,n)

                for (node_t j = child_ind_l; j < child_ind_u; ++j) {
                    if (elements[j].key < min_key) {
//                    if (elements[j].key < min_key || (elements[j].key >= min_key && elements[j].element < elements[i].element)) {
                        min_ind = j;
                        min_key = elements[j].key;
                    }
                }

                // Exchange?
                if (min_ind != i) {
                    swap(i, min_ind);
                    i = min_ind;
                } else {
                    break;
                }
            }
        }

        // Swap two elements in the heap.
        inline void swap(const node_t i, const node_t j) {
            element_t &el_i = elements[i];
            element_t &el_j = elements[j];

            // Exchange positions
            position[el_i.element] = j;
            position[el_j.element] = i;

            // Exchange elements
            element_t temp = el_i;
            el_i = el_j;
            el_j = temp;
        }

    private:

        // Number of elements in the heap.
        node_t n;

        // Number of maximal elements.
        node_t max_n;

        // Array of length heap_elements.
        std::vector<element_t> elements;

        // An array of positions for all elements.
        std::vector<node_t> position;
    };
}

