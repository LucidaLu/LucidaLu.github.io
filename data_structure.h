#ifndef LUCIDA_DATA_STRUCTURE
#define LUCIDA_DATA_STRUCTURE
#include "universal.h"
NAMESPACE_BEGIN(data_structure)
NAMESPACE_BEGIN(hashmap)
struct Node {
	int key,val;
	Node *pre;
	Node(int key,int val,Node *pre):key(key),val(val),pre(pre) {}
};
struct HashMap {
	int modu;
	Node **id;
	HashMap(int modu=99929):modu(modu),id((Node**)calloc(modu,sizeof(Node*))) {}
	~HashMap() {
		for(int i=0;i<modu;++i) {
			for(Node *p=id[i];p;) {
				Node *pre=p->pre;
				delete p;
				p=pre;
			}
		}
		free(id);
	}
	int Find(int key) {
		for(Node *p=id[key%modu];p;p=p->pre) if(p->key==key) {
			return p->val;
		}
		return -1;
	}
	void Insert(int key,int val) {
		int hs=key%modu;
		id[hs]=new Node(key,val,id[hs]);
	}
};
NAMESPACE_END
using hashmap::HashMap;

NAMESPACE_BEGIN(lct)
struct Node {
	Node *son[2];
	Data self,sum;
	Node() {}
};
NAMESPACE_END
using lct::LCT;
NAMESPACE_END
#endif
