#ifndef DSGVG_PHF_HPP
#define DSGVG_PHF_HPP

#include <cstdint>
#include "BooPHF.h"

namespace odgi {

typedef boomphf::SingleHashFunctor<uint64_t>  hasher_t;
typedef boomphf::mphf<uint64_t, hasher_t> boophf_t;

// iterator from disk file of uint64_t with buffered read,   todo template
class bfile_iterator : public std::iterator<std::forward_iterator_tag, uint64_t>{
	public:

	bfile_iterator()
	: _is(nullptr)
	, _pos(0) ,_inbuff (0), _cptread(0)
	{
		_buffsize = 10000;
		_buffer = (uint64_t *) malloc(_buffsize*sizeof(uint64_t));
	}

	bfile_iterator(const bfile_iterator& cr)
	{
		_buffsize = cr._buffsize;
		_pos = cr._pos;
		_is = cr._is;
		_buffer = (uint64_t *) malloc(_buffsize*sizeof(uint64_t));
		 memcpy(_buffer,cr._buffer,_buffsize*sizeof(uint64_t) );
		_inbuff = cr._inbuff;
		_cptread = cr._cptread;
		_elem = cr._elem;
	}

	bfile_iterator(FILE* is): _is(is) , _pos(0) ,_inbuff (0), _cptread(0)
	{
		_buffsize = 10000;
		_buffer = (uint64_t *) malloc(_buffsize*sizeof(uint64_t));
		int reso = fseek(_is,0,SEEK_SET);
		advance();
	}

	~bfile_iterator()
	{
		if(_buffer!=NULL)
			free(_buffer);
	}


	uint64_t const& operator*()  {  return _elem;  }

	bfile_iterator& operator++()
	{
		advance();
		return *this;
	}

	friend bool operator==(bfile_iterator const& lhs, bfile_iterator const& rhs)
	{
		if (!lhs._is || !rhs._is)  {  if (!lhs._is && !rhs._is) {  return true; } else {  return false;  } }
		assert(lhs._is == rhs._is);
		return rhs._pos == lhs._pos;
	}

	friend bool operator!=(bfile_iterator const& lhs, bfile_iterator const& rhs)  {  return !(lhs == rhs);  }
	private:
	void advance()
	{
		_pos++;

		if(_cptread >= _inbuff)
		{
			int res = fread(_buffer,sizeof(uint64_t),_buffsize,_is);
			_inbuff = res; _cptread = 0;

			if(res == 0)
			{
				_is = nullptr;
				_pos = 0;
				return;
			}
		}

		_elem = _buffer[_cptread];
		_cptread ++;
	}
	uint64_t _elem;
	FILE * _is;
	unsigned long _pos;

	uint64_t * _buffer; // for buffered read
	int _inbuff, _cptread;
	int _buffsize;
};


class file_binary{
	public:

	file_binary(const char* filename)
	{
		_is = fopen(filename, "rb");
		if (!_is) {
			throw std::invalid_argument("Error opening " + std::string(filename));
		}
	}

	~file_binary()
	{
		fclose(_is);
	}

	bfile_iterator begin() const
	{
		return bfile_iterator(_is);
	}

	bfile_iterator end() const {return bfile_iterator(); }

	size_t        size () const  {  return 0;  }//todo ?

	private:
	FILE * _is;
};

}

#endif
