#include <cmath>
#include <iostream>


// doubly linked list used to store items of a bucket
class item {

public: 

	int id;			// identifier
	double prio;	// priority 
	item *succ, *pred;
	 

	// constructor
	item(int _id, double _prio = 0) {
		this->id   = _id;
		this->prio = _prio;
		this->succ = NULL;
		this->pred = NULL;
	}

	// copy constructor 
	item(const item &_obj) {
		this->id =   _obj.id;
		this->prio = _obj.prio;
		this->succ = _obj.succ;
		this->pred = _obj.pred;		
	}

	// desctructor
	~item() { 
		item *curr = this;

		while (curr != NULL) {
			item *tmp = curr;
			curr = curr->succ;
			free(tmp);
		}    	
	}

    void print() {
		item *curr = this;

		while (curr != NULL) {
			printf("%d:%.2f ", curr->id, curr->prio);
			curr = curr->succ;
		}    	
		std::cout << std::flush;
    }

};


class bucket_list {

public: 

	int    max_items;   // maximum number of items
	int    size; 		// current number of items in bucket lists

	int    cap; 		// number of buckets (fixed in constructor)

	double base; 	    // original base value (as set initially)
	double beta; 	    // scalar

	double cur_base;	// current base value (updated throughout)
	int    idx_0;		// largest index of unused bucket (0 initially)
						// invariant: buckets 0, ..., idx_0 all empty
						// relation: cur_base = beta^idx_0 base

	item   **bucket; 	// array of pointers to head of bucket lists
	item   **map;		// array of pointers to items of id
	
	bool is_init;


	bucket_list() {
		// initialize to some temporary values
		this->is_init = false;
		this->max_items = 0;
		this->size  = 0;
		this->cap  = 0; 
		this->base = 0.0;
		this->beta = 0.0; 
		this->cur_base = 0.0;
		this->idx_0 = 0;
		this->bucket = NULL;
		this->map = NULL;
	}
	
	void set_base(double _base){
		this->base = _base;
	}
	
	void init(int _max_items, double _base, double _beta, int _cap = 1){
		this->is_init = true;
		this->max_items = _max_items;
		this->size  = 0;
		this->cap  = _cap; 
		this->base = _base;
		this->beta = _beta; 
		this->cur_base = _base;		
		this->idx_0 = 0;
		this->bucket  = new item* [_cap];
		this->map = new item* [_max_items];			

		// initialize pointers
		for (int i = 0; i < _cap; i++) this->bucket[i] = NULL;				
		for (int i = 0; i < _max_items; i++) this->map[i] = NULL;			
	}


	void stats() {
		printf("\n----- %d buckets, %d size, %d idx_0 -----\n", cap, size, idx_0);
		for(int i = idx_0+1; i < cap; i++) {
			printf("\nbucket %d: ", i);
			this->bucket[i]->print();
		}	
		printf("\n\n");
	}

	int get_bucket_number(double _prio) {
		// bucket 0: 				prio <= base
		// bucket 1: 		base < 	prio <= beta   * base
		// bucket 2: beta * base < 	prio <= beta^2 * base, ...

		int idx = (_prio > base ? ceil( log( _prio / base ) / log( beta ) ) : 0);

		if (idx > cap) {
			std::cout << "\n error: bucket_list: priority out of range\n\n";
			return 0;
		}
		return idx;
	}


	void insert(int _id, double _prio) {
		// insert item _id at front of bucket asociated with _prio
		int idx = get_bucket_number(_prio); 

		if (idx <= idx_0) return;  // item irrelevant 
	
		// create a new item
		size++;
		item *new_item = new item(_id, _prio);
		this->map[_id] = new_item;

		// add new_item to corresonding bucket
		if (this->bucket[idx] == NULL) {
			// create new bucket list 
			this->bucket[idx] = new_item;
		}
		else {
			// add new element at front of existing bucket list
			item *head = this->bucket[idx];
			new_item->succ = head;
			head->pred = new_item; 
			this->bucket[idx] = new_item;
		}
	}

	
	void move(int _id, double _new_prio) {
		// moves item _id to bucket associated with new priority _new_prio

		if (this->map[_id] == NULL) return;

		item *cur_item = this->map[_id];
		double cur_prio = cur_item->prio;
		int cur_idx = get_bucket_number(cur_prio);
		int new_idx = get_bucket_number(_new_prio); 
		
		if ( new_idx <= idx_0 || (cur_idx == new_idx) ) return;  // bucket irrelevant or remains unchanged

		// get pred and succ item of cur_item
		item *tmp_pre = cur_item->pred; 
		item *tmp_suc = cur_item->succ; 

		// redirect pointers of neighboring items
		if (tmp_pre != NULL) tmp_pre->succ = tmp_suc;
		if (tmp_suc != NULL) tmp_suc->pred = tmp_pre;

		// check if head of bucket referred to cur_item
		if (this->bucket[cur_idx] == cur_item) this->bucket[cur_idx] = tmp_suc;

		// update info of cur_item
		cur_item->prio = _new_prio;
		cur_item->succ = NULL;
		cur_item->pred = NULL;

		// insert cur_item at front of bucket[new_idx]
		if (this->bucket[new_idx] == NULL) {
			// create new bucket list 
			this->bucket[new_idx] = cur_item;
		}
		else {
			// add new element at front of existing bucket list
			item *head = this->bucket[new_idx];
			cur_item->succ = head;
			head->pred = cur_item; 
			this->bucket[new_idx] = cur_item;
		}
	}


	// increases idx_0 until first non-empty bucket is found
	int advance(double q_front) { 

		int tmp = idx_0;
		while (empty(idx_0) and (q_front > (cur_base)) ) {
			cur_base *= beta;
			idx_0++;
			if (idx_0 == cap) return 0; // reached end of buckets
		} 
		return idx_0 - tmp; // return number of update steps
	}

	bool contains(int _id) { return ((this->is_init) and (this->map[_id] != NULL)); }

	// checks whether all buckets are empty 
	bool empty() { return (size == 0); } 

	// checks whether bucket[_idx] is empty
	bool empty(int _idx) { 
		return (this->bucket[_idx] == NULL); 
	} 

	int front(int _idx) {
		// returns the id of the head element of bucket[_idx] 
		if (empty(_idx)) return 0; 
		else return this->bucket[_idx]->id;
	}

	int pop(int _idx) {
		// removes first item from bucket[_idx] and returns its id

		if (empty(_idx)) return 0; 

		item *cur_item = this->bucket[_idx];
		int cur_id = cur_item->id;	
		item *tmp_succ = cur_item->succ; 

		// redirect pred pointer of tmp_succ to NULL (head of list)
		if (tmp_succ) tmp_succ->pred = NULL;

		// set head of bucket to new head item (or NULL if empty)
		this->bucket[_idx] = tmp_succ;

		// remove current item
		this->map[cur_id] = NULL;
		free(cur_item);		
		size--;

		return cur_id;
	}


	// removes first item from bucket[idx_0] and returns its id
	int pop_0() { return pop(idx_0); }


	void remove(int _id) {
		// removes item _id from respective bucket list

		if (this->map[_id] == NULL) return;

		item *cur_item = this->map[_id];
		double cur_prio = cur_item->prio;	
		int cur_idx = get_bucket_number(cur_prio);

		// get pred and succ item of cur_item
		item *tmp_pre = cur_item->pred; 
		item *tmp_suc = cur_item->succ; 

		// redirect pointers of neighboring items
		if (tmp_pre != NULL) tmp_pre->succ = tmp_suc;
		if (tmp_suc != NULL) tmp_suc->pred = tmp_pre;

		// check if head of bucket referred to cur_item
		if (this->bucket[cur_idx] == cur_item) this->bucket[cur_idx] = tmp_suc;

		// remove cur_item 
		this->map[_id] = NULL;
		free(cur_item);		
		size--;
	}




};