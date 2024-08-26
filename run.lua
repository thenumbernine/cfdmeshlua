#!/usr/bin/env luajit
require 'ext'
local gl = require 'gl'
local ig = require 'imgui'
local matrix = require 'matrix'

local App = require 'imguiapp.withorbit'()
App.viewUseGLMatrixMode = true
App.title = 'cfd mesh'

local Vertex = class()
function Vertex:init(x,y)
	self.pos = matrix{x,y}
	self.edges = table()
end

local Edge = class()
function Edge:init(...)
	self.vtxs = table{...}
	self.cells = table()
end

local Cell = class()
function Cell:init()
	self.edges = table()
	self.vtxs = table()
end

local vtxs = table()
local edges = table()
local cells = table()

local mergeUniqueVtxs = false

local function addVtx(x)
	if mergeUniqueVtxs then
		for _,v in ipairs(vtxs) do
			if v.pos == x then return v end
		end
	end
	local v = Vertex()
	v.pos = matrix(x)
	vtxs:insert(v)
	return v
end

local function addEdge(i,j)
	local a,b = vtxs[i], vtxs[j]
	for _,e in ipairs(edges) do
		if (e.vtxs[1] == a and e.vtxs[2] == b)
		or (e.vtxs[1] == b and e.vtxs[2] == a)
		then
			return e
		end
	end
	local e = Edge(a,b)
	edges:insert(e)
	a.edges:insert(e)
	b.edges:insert(e)
	return e
end

local function addCell(...)
	local c = Cell()
	local n = select('#', ...)
	local indexes = table()
	for i=1,n do
		local vtxindex = select(i, ...)
		local v = assert(vtxs[vtxindex], "couldn't find vertex "..vtxindex)
		c.vtxs:insert(v)
		indexes:insert(vtxindex)
	end
	for i=1,n do
		c.edges:insert(addEdge(indexes[i], indexes[i%n+1]))
	end
	cells:insert(c)
	c.pos = c.vtxs:map(function(v) return v.pos end):sum() / n
	for _,e in ipairs(c.edges) do
		e.cells:insert(c)
	end
	return c
end

local function polyVol(...)
	local n = select('#', ...)
	local v = 0
	for i=1,n do
		local pi = select(i, ...)
		local pj = select(i%n+1, ...)
		v = v + .5 * (pi[1] * pj[2] - pi[2] * pj[1])
	end
	return v
end

local function calc_hTotal(rho, P, ETotal)
	return (ETotal + P) / rho
end

local heatCapacityRatio = 1.4
local function calcSpeedOfSound(vx, vy, vz, hTotal)
	return math.sqrt((heatCapacityRatio - 1) * (hTotal - .5 * (vx*vx + vy*vy + vz*vz)))
end

local function consFromPrim(W)
	local rho, vx, vy, vz, P = W:unpack()
	return matrix{
		rho,
		rho * vx,
		rho * vy,
		rho * vz,
		P/(heatCapacityRatio-1) + .5*rho*(vx*vx + vy*vy + vz*vz),
	}
end

local function primFromCons(U)
	local rho, mx, my, mz, ETotal = U:unpack()
	local vx = mx / rho
	local vy = my / rho
	local vz = mz / rho
	local P = (heatCapacityRatio - 1) * (ETotal - .5 * rho * (vx*vx + vy*vy + vz*vz))
	return matrix{rho, vx, vy, vz, P}
end

local mins = matrix{-1,-1}
local maxs = matrix{1,1}
function App:initGL(...)
	App.super.initGL(self, ...)
	self.view.ortho = true
	
	-- TODO base these on the grid bounds
	self.view.orthoSize = 1.5
	self.view.znear = 1
	self.view.zfar = 1000

--[[ grid
	local m = 30
	local n = 20
	for i=1,m+1 do
		for j=1,n+1 do
			addVtx(matrix{
				(i-.5)/(m+1)*(maxs[1] - mins[1]) + mins[1],
				(j-.5)/(n+1)*(maxs[2] - mins[2]) + mins[2]})
		end
	end

	for i=1,m do
		for j=1,n do
			local c = addCell(
				1 + j-1 + (n+1) * (i-1),
				1 + j-1 + (n+1) * i,
				1 + j + (n+1) * i,
				1 + j + (n+1) * (i-1))
			if c.pos[1] < 0 and c.pos[2] < 0 then
				c.U = consFromPrim(matrix{
					1,
					0, 0, 0,
					1
				})
			else
				c.U = consFromPrim(matrix{
					.125,
					0,0,0,
					.1
				})
			end
		end
	end
--]]


-- [[ file
	local ls = assert(path'grids/n0012_113-33.p2dfmt':read()):trim():split'\n'
	local first = ls:remove(1)
	local m, n = ls:remove(1):trim():split'%s+':map(function(l) return tonumber(l) end):unpack()
	local x = ls:concat():trim():split'%s+':map(function(l) return tonumber(l) end)
	assert(#x == 2*m*n)
	print(m, n, m*n)
	-- [[
	local us = x:sub(1,m*n)
	local vs = x:sub(m*n+1)
	assert(#us == #vs)
	for i=1,#us do
		local u,v = us[i], vs[i]
		print(u,v)
		addVtx(matrix{u,v})
	end
	--]]
	--[[
	local k = 1
	for i=1,m*n do
		local u = x[k] k=k+1
		local v = x[k] k=k+1
		print(u,v)
		addVtx(matrix{u,v})
	end
	--]]
	assert(#vtxs == m*n, "expected "..#vtxs.." to equal "..(m*n))
	
	for i=1,n-1 do
		for j=1,m-1 do
			local c = addCell(
				1 + j-1 + m * (i-1),
				1 + j-1 + m * i,
				1 + j + m * i,
				1 + j + m * (i-1))
			
			c.U = consFromPrim(matrix{
				1,
				.1, 0, 0,
				1
			})
		end
	end
--]]


	for _,e in ipairs(edges) do
		local a,b = e.vtxs:unpack()
		e.pos = (a.pos + b.pos) * .5
		e.delta = a.pos - b.pos
		e.length = e.delta:norm()
		e.normal = matrix{-e.delta[2], e.delta[1]}
		e.normal = e.normal / e.normal:norm()
		
		local a,b = e.cells:unpack()
		if a and b then
			-- make sure the normal points to a
			if (a.pos - e.pos):dot(e.normal) < 0 then
				a,b = b,a
				e.cells = table{a,b}
			end
			e.cellDist = (b.pos - a.pos):norm()
		else
			local c = a or b
			-- for ghost state's sake:
			e.cellDist = (c.pos - e.pos):norm() * 2
		end
	end

	for _,c in ipairs(cells) do
		c.volume = polyVol(c.vtxs[1].pos, c.vtxs[2].pos, c.vtxs[3].pos, c.vtxs[4].pos)
	end
end

-- rotate vx,vy such that n now points along the x dir
local function rotateTo(vx, vy, n)
	return vx * n[1] + vy * n[2], vy * n[1] - vx * n[2]
end

-- rotate vx,vy such that the x dir now points along n 
local function rotateFrom(vx, vy, n)
	return vx * n[1] - vy * n[2], vy * n[1] + vx * n[2]
end

running = false
showVtxs = true
showCellCenters = false
showEdges = true

function App:draw()
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	gl.glBegin(gl.GL_QUADS)
	for _,c in ipairs(cells) do
		local W = primFromCons(c.U)
		gl.glColor3f(
			W[1],
			math.sqrt(W[2]*W[2] + W[3]*W[3] + W[4]*W[4]),
			W[5]
		)
		for _,v in ipairs(c.vtxs) do
			gl.glVertex2d(v.pos:unpack())
		end
	end
	gl.glEnd()

	if showVtxs or showCellCenters then
		gl.glColor3f(1,1,1)
		gl.glPointSize(3)
		gl.glBegin(gl.GL_POINTS)
		if showVtxs then
			for _,v in ipairs(vtxs) do
				gl.glVertex2d(v.pos:unpack())
			end
		end
		if showCellCenters then
			for _,c in ipairs(cells) do
				gl.glVertex2d(c.pos:unpack())
			end
		end
		gl.glEnd()
		gl.glPointSize(1)
	end
	if showEdges then
		gl.glBegin(gl.GL_LINES)
		for _,e in ipairs(edges) do
			local a,b = e.vtxs:unpack()
			gl.glVertex2d(a.pos:unpack())
			gl.glVertex2d(b.pos:unpack())
		end
		gl.glEnd()
	end
end

function App:step()
	-- calculate dt
	local result = math.huge
	local cfl = .5
	for _,c in ipairs(cells) do
		local W = primFromCons(c.U)
		local rho,vx,vy,vz,P = W:unpack()
		local Cs = math.sqrt(heatCapacityRatio * P / rho)
		for j=2,4 do
			local v = W[j]
			--local vx, vy = rotateTo(vx,vy, e.normal)
			--assert(vy == 0)
			local lambdaMin = v - Cs
			local lambdaMax = v + Cs
			lambdaMin = math.min(0, lambdaMin)
			lambdaMax = math.max(0, lambdaMax)
			-- TODO a better way to do this.  maybe use edges' lambdas?  maybe do this after calculating the eigenbasis?
			local dx = math.sqrt(c.volume)
			local dum = dx / (math.abs(lambdaMax - lambdaMin) + 1e-9)
			result = math.min(result, dum)
		end
	end
	local dt = result * cfl

	local function getEdgeStates(e)
		local cL, cR = e.cells:unpack()
		local UL, UR
		if cL and cR then
			UL = matrix(assert(cL.U))
			UR = matrix(assert(cR.U))
		elseif cL then
			UL = matrix(assert(cL.U))
			UR = matrix(UL)
			local m = matrix{UR[2],UR[3]}
			m = m - e.normal * (2 * (e.normal * m))
			UR[2], UR[3] = m:unpack()
		elseif cR then
			UR = matrix(assert(cR.U))
			UL = matrix(UR)
			local m = matrix{UL[2],UL[3]}
			m = m - e.normal * (2 * (e.normal * m))
			UL[2], UL[3] = m:unpack()
		else
			error"here"
		end
		return UL, UR
	end

	for _,e in ipairs(edges) do
		local UL, UR = getEdgeStates(e)

			-- 1) roe values at edge 
		do
			-- rotate to align edge normal to x axis
			-- so x-direction flux jacobian is good for calculating the flux 
			UL[2], UL[3] = rotateTo(UL[2], UL[3], e.normal)
			UR[2], UR[3] = rotateTo(UR[2], UR[3], e.normal)

			local ETotalL = UL[5]
			local rhoL, vxL, vyL, vzL, PL = primFromCons(UL):unpack()
			local hTotalL = calc_hTotal(rhoL, PL, ETotalL)
			local sqrtRhoL = math.sqrt(rhoL)
			
			local ETotalR = UR[5]
			local rhoR, vxR, vyR, vzR, PR = primFromCons(UR):unpack()
			local hTotalR = calc_hTotal(rhoR, PR, ETotalR)
			local sqrtRhoR = math.sqrt(rhoR)
			
			local rho = sqrtRhoL * sqrtRhoR
			local vx = (sqrtRhoL * vxL + sqrtRhoR * vxR) / (sqrtRhoL + sqrtRhoR)
			local vy = (sqrtRhoL * vyL + sqrtRhoR * vyR) / (sqrtRhoL + sqrtRhoR)
			local vz = (sqrtRhoL * vzL + sqrtRhoR * vzR) / (sqrtRhoL + sqrtRhoR)
			local hTotal = (sqrtRhoL * hTotalL + sqrtRhoR * hTotalR) / (sqrtRhoL + sqrtRhoR)
			local Cs = calcSpeedOfSound(vx, vy, vz, hTotal)
			local CsSq = Cs * Cs
			--e.roe = matrix{rho, vx, vy, vz, hTotal, Cs}
		
			-- 2) eigenbasis at interface
		
			local vSq = vx*vx + vy*vy + vz*vz

			local lambdas = matrix{vx - Cs, vx, vx, vx, vx + Cs}
			local gamma_1 = heatCapacityRatio - 1	
			local evL = matrix{
				{(.5 * gamma_1 * vSq + Cs * vx) / (2 * CsSq),	(-Cs - gamma_1 * vx) / (2 * CsSq),	-gamma_1 * vy / (2 * CsSq),		-gamma_1 * vz / (2 * CsSq),	gamma_1 / (2 * CsSq),	},
				{1 - gamma_1 * vSq / (2 * CsSq),				gamma_1 * vx / CsSq,				gamma_1 * vy / CsSq,			gamma_1 * vz / CsSq,		-gamma_1 / CsSq,		},
				{-vy,											0,									1,								0,							0,						}, 
				{-vz,											0,									0,								1,							0,						},
				{(.5 * gamma_1 * vSq - Cs * vx) / (2 * CsSq),	(Cs - gamma_1 * vx) / (2 * CsSq),	-gamma_1 * vy / (2 * CsSq),		-gamma_1 * vz / (2 * CsSq),	gamma_1 / (2 * CsSq),	},
			}
	
			local evR = matrix{
				{1, 				1, 			0,		0,		1,				},
				{vx - Cs, 			vx, 		0,		0,		vx + Cs,		},
				{vy,				vy,			1,		0,		vy,				},
				{vz,				vz,			0,		1,		vz,				},
				{hTotal - Cs * vx, .5 * vSq, 	vy,		vz,		hTotal + Cs * vx},
			}
		
			local dU = UR - UL
			local dUTilde = evL * dU
		
			local fluxTilde = matrix()
			for j=1,5 do
				local lambda = lambdas[j]
				local phi = 0
				local sgnLambda = lambda >= 0 and 1 or -1
				local dx = e.cellDist
				local epsilon = lambda * dt / dx
				fluxTilde[j] = -.5 * lambda * dUTilde[j] * (sgnLambda + phi * (epsilon - sgnLambda))
			end
		
			local UAvg = (UR + UL) * .5
			local UAvgTilde = evL * UAvg
			fluxTilde = fluxTilde + lambdas:emul(UAvgTilde)
		
			local flux = evR * fluxTilde
			-- here's the flux, aligned along the normal
		
			flux[2], flux[3] = rotateFrom(flux[2], flux[3], e.normal)
			-- here's the flux in underlying coordinates
			e.flux = flux
		end
	end

	for _,c in ipairs(cells) do
		local dU_dt = matrix{0,0,0,0,0}
		for _,e in ipairs(c.edges) do
			if e.flux then
				if c == e.cells[1] then
					dU_dt = dU_dt - e.flux * (e.length / c.volume)
				else
					dU_dt = dU_dt + e.flux * (e.length / c.volume)
				end
			end
		end
		c.U = c.U + dU_dt * dt
	end
end

function App:update()
	self:draw()
	if running then
		self:step()
	end
	App.super.update(self)
end

function App:updateGUI()
	ig.luatableCheckbox('running', _G, 'running')
	ig.luatableCheckbox('showVtxs', _G, 'showVtxs')
	ig.luatableCheckbox('showCellCenters', _G, 'showCellCenters')
	ig.luatableCheckbox('showEdges', _G, 'showEdges')
	ig.luatableCheckbox('ortho', self.view, 'ortho')
end

return App():run()
