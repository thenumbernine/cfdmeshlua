#!/usr/bin/env luajit
require 'ext'
local gl = require 'gl'
local ImGuiApp = require 'imguiapp'
local View = require 'glapp.view'
local Orbit = require 'glapp.orbit'
local App = class(Orbit(View.apply(ImGuiApp)))
local matrix = require 'matrix'
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

local function addVtx(x)
	for _,v in ipairs(vtxs) do
		if v.pos == x then return v end
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

local function addCell(a,b,c,d)
	local v = Cell()
	v.vtxs:insert(vtxs[a])
	v.vtxs:insert(vtxs[b])
	v.vtxs:insert(vtxs[c])
	v.vtxs:insert(vtxs[d])
	v.edges:insert(addEdge(a,b))
	v.edges:insert(addEdge(b,c))
	v.edges:insert(addEdge(c,d))
	v.edges:insert(addEdge(d,a))
	cells:insert(v)
	v.pos = v.vtxs:map(function(v) return v.pos end):sum() * .25
	for _,e in ipairs(v.edges) do
		e.cells:insert(v)
	end
	return v
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
function App:initGL()
	self.view.ortho = true
	self.view.orthoSize = 1.5
	local n = 50
	for i=1,n+1 do
		for j=1,n+1 do
			addVtx(matrix{
				(i-.5)/(n+1)*(maxs[1] - mins[1]) + mins[1],
				(j-.5)/(n+1)*(maxs[2] - mins[2]) + mins[2]})
		end
	end
	for i=1,n do
		for j=1,n do
			local c = addCell(
				1 + i-1 + (n+1) * (j-1),
				1 + i-1 + (n+1) * j,
				1 + i + (n+1) * j,
				1 + i + (n+1) * (j-1))
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

function App:update()
	gl.glClear(gl.GL_COLOR_BUFFER_BIT)
	App.super.update(self)
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
	gl.glColor3f(1,1,1)
	gl.glPointSize(3)
	gl.glBegin(gl.GL_POINTS)
	for _,v in ipairs(vtxs) do
		gl.glVertex2d(v.pos:unpack())
	end
	for _,c in ipairs(cells) do
		gl.glVertex2d(c.pos:unpack())
	end
	gl.glEnd()
	gl.glPointSize(1)
	gl.glBegin(gl.GL_LINES)
	for _,e in ipairs(edges) do
		local a,b = e.vtxs:unpack()
		gl.glVertex2d(a.pos:unpack())
		gl.glVertex2d(b.pos:unpack())
	end
	gl.glEnd()

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

	for _,e in ipairs(edges) do
		local cL, cR = e.cells:unpack()
		if cL and cR then

			-- 1) roe values at edge 
			
			local UL = matrix(assert(cL.U))
			local UR = matrix(assert(cR.U))

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

App():run()
