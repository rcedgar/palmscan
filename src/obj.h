#ifndef obj_h
#define obj_h

#define TRACE_OBJS			0
#define TRACK_OBJ_THREAD	0

#include "objtype.h"

/***
Objects are created and deleted only by ObjMgr.

Objects can contain other objects.

Getting new object:
	ptr = ObjMgr::GetXxx()
	New object has ref count = 1.

Adding reference:
	Convention: one variable per counted reference.
    assert(var == 0)
	var = ptr
	OM->Up(var)
	Up() is similar to IUnknown::AddRef().

Releasing reference:
	OM->Down(var)
	var = 0
	Down() is similar to IUknown::Release().

When objects contain references to other objects:
	In Init/Create()-like function, m_xxx = OM->GetXxx().
	In OnZeroRefCount(), call OM->Down(m_xxx).

There is one ObjMgr per thread, avoids need for mutexes etc.
No Clear() member, this is confusing.

Memory allocation convention (not enforced or expressed):
	Objects use grow-only memory strategy.
	D'tor frees memory.
    D'tor is the only place memory is freed (except growing).
***/

class ObjMgr;

class Obj
	{
	friend class ObjMgr;
#if	TRACK_OBJ_THREAD
public:
	int m_OMPThreadIndex;
#endif

protected:
	ObjType m_Type;
	ObjMgr *m_ObjMgr;
	unsigned m_RefCount;
	Obj *m_Fwd;
	Obj *m_Bwd;
#if	TRACE_OBJS
	unsigned m_UID;
	const char *m_SourceFileName;
	unsigned m_SourceLineNr;
#endif

private:
	Obj();

protected:
	Obj(ObjType Type, ObjMgr *OM)
		{
		m_Type = Type;
		m_ObjMgr = OM;
		m_RefCount = 0;
		m_Fwd = 0;
		m_Bwd = 0;
		}

	virtual ~Obj()
		{
		}

public:
	unsigned GetRefCount() const
		{
		return m_RefCount;
		}

	ObjMgr *GetObjMgr()
		{
		return m_ObjMgr;
		}

// Override if contains objects, must Down() them.
	virtual void OnZeroRefCount() {}
	virtual unsigned GetMemBytes() const = 0;
	};

#if	TRACE_OBJ
#define	Up()	Up_(__FILE__, __LINE__)
#define Down()	Down_(__FILE__, __LINE__)
#endif

#endif // obj_h
