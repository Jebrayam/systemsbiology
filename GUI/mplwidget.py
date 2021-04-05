# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 16:49:03 2020

@author: Usuario
"""

#from PyQt5.QtWidgets import*
from PyQt5.QtWidgets import QWidget, QVBoxLayout
from matplotlib.backends.backend_qt5agg import FigureCanvas
from matplotlib.figure import Figure

    
class MplWidget(QWidget):
    
    def __init__(self, parent = None):

        QWidget.__init__(self, parent)
        
        self.canvas = FigureCanvas(Figure())
        
        vertical_layout = QVBoxLayout()
        vertical_layout.addWidget(self.canvas)
        
        self.canvas.axes = self.canvas.figure.add_subplot(111)
        self.canvas.axes.tick_params(axis='both', direction='in', labelsize=8)
        self.setLayout(vertical_layout)