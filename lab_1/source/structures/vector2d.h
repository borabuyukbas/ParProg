#pragma once

#include <exception>

template <typename _Ty>
class Vector2d
{
public:

	Vector2d() :
		m_x(0),
		m_y(0)
	{}

	Vector2d(_Ty p_x, _Ty p_y) :
		m_x(p_x), m_y(p_y)
	{}

	_Ty operator[](const int& index)
	{
		// exception handling
		if (index < 0 || index > 1)
			throw(std::exception{});

		return index == 0 ? m_x : m_y;
	}


	Vector2d<_Ty> operator+(Vector2d vec)
	{
		return Vector2d<_Ty>(this->m_x + vec.m_x, this->m_y + vec.m_y);
	}


	Vector2d<_Ty> operator-(Vector2d vec)
	{
		return Vector2d<_Ty>(this->m_x - vec.m_x, this->m_y - vec.m_y);
	}


	Vector2d<_Ty> operator*(Vector2d vec)
	{
		return Vector2d<_Ty>(this->m_x * vec.m_x, this->m_y * vec.m_y);
	}


	Vector2d<_Ty> operator/(Vector2d vec)
	{
		return Vector2d<_Ty>(this->m_x / vec.m_x, this->m_y / vec.m_y);
	}


	bool operator==(Vector2d vec)
	{
		return ((this->m_x == vec.m_x) && (this->m_y == vec.m_y));
	}


	void set(_Ty p_x, _Ty p_y)
	{
		this[0] = p_x;
		this[1] = p_y;
	}

private:
	_Ty m_x, m_y;

};



