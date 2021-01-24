//
// Created by af on 17.01.21.
//

#include "aml_lua_binding.h"
#include <iostream>

int amlAdd4xN(lua_State *L) { // TODO detect the length of the table
	if (lua_istable(L, 1) && lua_istable(L, 2)) {

		double a0 = 0;
		double a1 = 0;
		double a2 = 0;
		double a3 = 0;
		double b0 = 0;
		double b1 = 0;
		double b2 = 0;
		double b3 = 0;
		lua_pushinteger(L, 1);
		int type = lua_gettable(L, -2);
		if (type == LUA_TNUMBER) {
			a0 = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);

		//

		lua_pushinteger(L, 2);
		type = lua_gettable(L, -2);
		if (type == LUA_TNUMBER) {
			a1 = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);

		//

		lua_pushinteger(L, 3);
		type = lua_gettable(L, -2);
		if (type == LUA_TNUMBER) {
			a2 = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);

		//
		lua_pushinteger(L, 4);
		type = lua_gettable(L, -2);
		if (type == LUA_TNUMBER) {
			a3 = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);

		lua_pop(L, 1);
		//
		//

		lua_pushinteger(L, 1);
		type = lua_gettable(L, -2);
		if (type == LUA_TNUMBER) {
			b0 = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);

		//

		lua_pushinteger(L, 2);
		type = lua_gettable(L, -2);
		if (type == LUA_TNUMBER) {
			b1 = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);

		//

		lua_pushinteger(L, 3);
		type = lua_gettable(L, -2);
		if (type == LUA_TNUMBER) {
			b2 = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);

		//
		lua_pushinteger(L, 4);
		type = lua_gettable(L, -2);
		if (type == LUA_TNUMBER) {
			b3 = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
		lua_pop(L, 1);

		//


		VectorDouble4D a(b0, b1, b2, b3);//I know confusing
		VectorDouble4D b(a0, a1, a2, a3);//doesn't mater for addition but subtraction

		a += b;

		a0 = a[0];
		a1 = a[1];
		a2 = a[2];
		a3 = a[3];

		lua_newtable(L);
		lua_pushnumber(L, a0);
		lua_seti(L, -2, 1);

		lua_pushnumber(L, a1);
		lua_seti(L, -2, 2);
		lua_pushnumber(L, a2);
		lua_seti(L, -2, 3);
		lua_pushnumber(L, a3);
		lua_seti(L, -2, 4);

		return 1;
	} else {
		std::cout << "no table given" << std::endl;
		return 0;
	}
}

int amlMap4xN(lua_State *L) {
	if (lua_istable(L, -1) && lua_isnumber(L, -2) && lua_isnumber(L, -3) && lua_isnumber(L, -4) &&
		lua_isnumber(L, -5)) {

		double a0 = 0;
		double a1 = 0;
		double a2 = 0;
		double a3 = 0;
		lua_pushinteger(L, 1);
		int type = lua_gettable(L, -2);
		if (type == LUA_TNUMBER) {
			a0 = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);

		//

		lua_pushinteger(L, 2);
		type = lua_gettable(L, -2);
		if (type == LUA_TNUMBER) {
			a1 = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);

		//

		lua_pushinteger(L, 3);
		type = lua_gettable(L, -2);
		if (type == LUA_TNUMBER) {
			a2 = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);

		//
		lua_pushinteger(L, 4);
		type = lua_gettable(L, -2);
		if (type == LUA_TNUMBER) {
			a3 = lua_tonumber(L, -1);
		}
		lua_pop(L, 1);
		lua_pop(L, 1);
		//
		//
		VectorDouble4D a(a0, a1, a2, a3);

		double lowerInput = lua_tonumber(L, -3);
		double upperInput = lua_tonumber(L, -4);
		double lowerOutput = lua_tonumber(L, -1);
		double upperOutput = lua_tonumber(L, -2);

		a.map(upperInput, lowerInput, upperOutput, lowerOutput);

		a0 = a[0];
		a1 = a[1];
		a2 = a[2];
		a3 = a[3];

		lua_newtable(L);
		lua_pushnumber(L, a0);
		lua_seti(L, -2, 1);

		lua_pushnumber(L, a1);
		lua_seti(L, -2, 2);
		lua_pushnumber(L, a2);
		lua_seti(L, -2, 3);
		lua_pushnumber(L, a3);
		lua_seti(L, -2, 4);

		return 1;
	} else {
		std::cout << "no table given" << std::endl;
		return 0;
	}
}


void initAmlLua(lua_State *L) {
	lua_register(L, "addVec4", amlAdd4xN);
	lua_register(L, "mapVec4", amlMap4xN);
}