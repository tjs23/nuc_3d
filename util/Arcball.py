from math import acos, sqrt, pi
import numpy
from numpy.linalg import norm

from . import XForm

class ArcBall:
	def __init__(self):
		self.old_rotation = numpy.identity(4)
		self.points = [None, None]
	
	def surfacePoint(self, *pt):
		pt = numpy.array(pt)

		surface_pt = numpy.zeros(3)
		surface_pt[0:2] = pt
		if norm(pt) < 1:
			surface_pt[2] = sqrt(1 - (pt ** 2).sum())
		else:
			surface_pt[0:2] = pt / norm(pt) 

		return surface_pt
	
	def startRotation(self, *pt):
		self.points[0] = self.surfacePoint(*pt)

	def updateRotation(self, *pt):
		self.points[1] = self.surfacePoint(*pt)

	def finishRotation(self):
		self.old_rotation = self.totalRotation()
		self.points = [None, None]
		
	def currentRotation(self):
		if None in self.points:
			theta = 0
		else:
			# FIXME: ValueError: math domain error here. Not sure what values were in self.points.
			theta = acos(numpy.dot(*self.points))

		if theta == 0:
			# n is arbitrary as there is no rotation
			n = (1, 0, 0)
		elif theta == pi:
			# self.points are on a line, need to be orthogonal to either
			# could do branchless as input is guaranteed to be normalized.
			p = self.points[0]
			if abs(p[0]) > abs(p[1]):
				n = (-p[1], p[0], 0)
			else:
				n = (0, -p[2], p[1])
		else:
			n = numpy.cross(*self.points)

		return XForm.rotateAngleAxis(theta, n)
	
	def totalRotation(self):
		return self.currentRotation().dot(self.old_rotation)
