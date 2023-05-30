"""
small widget to control wind settings presets taken from VehiclePhysicalConstants
"""

from PyQt5.QtCore import *
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
import ece163.Constants.VehiclePhysicalConstants as VehiclePhysicalConstants

widgetName = "Wind Control"

class WindControl(QWidget):
	def __init__(self, aeroInstance, parent=None):
		"""
		generates a wind control widget using the lists in phsyical parameters. Needs a simulate instance so it can set them.

		:param simulateInstance: simulate instance used for setting the wind values
		"""

		super().__init__(parent)

		self.aeroInstance = aeroInstance

		usedLayout = QVBoxLayout()
		self.setLayout(usedLayout)

		steadyWindwsColumn = QVBoxLayout()
		gustWindwsColumn = QVBoxLayout()

		windSelectionBox = QHBoxLayout()
		usedLayout.addLayout(windSelectionBox)
		windSelectionBox.addLayout(steadyWindwsColumn)
		windSelectionBox.addLayout(gustWindwsColumn)

		steadyWindwsColumn.addWidget(QLabel("Steady Windws"))
		gustWindwsColumn.addWidget(QLabel("Gust Winds"))

		self.steadyWinds = list()
		self.steadyButtonsGroup = QButtonGroup()
		for index, (name, value) in enumerate(VehiclePhysicalConstants.SteadyWinds):
			newRadio = QRadioButton(name)
			self.steadyWinds.append(value)
			steadyWindwsColumn.addWidget(newRadio)
			self.steadyButtonsGroup.addButton(newRadio, index)

		self.steadyButtonsGroup.button(0).setChecked(True)
		# print(self.steadyButtonsGroup.checkedId())

		self.gustWinds = list()
		self.gustButtonsGroup = QButtonGroup()
		for index, (name, value) in enumerate(VehiclePhysicalConstants.GustWinds):
			newRadio = QRadioButton(name)
			self.gustWinds.append(value)
			gustWindwsColumn.addWidget(newRadio)
			self.gustButtonsGroup.addButton(newRadio, index)

		self.gustButtonsGroup.button(0).setChecked(True)

		applyWindsButton = QPushButton("Apply Winds")
		applyWindsButton.clicked.connect(self.applyWindsResponse)
		usedLayout.addWidget(applyWindsButton)

		steadyWindwsColumn.addStretch()
		gustWindwsColumn.addStretch()
		return

	def applyWindsResponse(self):
		"""
		Reads in current wind radio buttons and adjusts the underlying model to the new parameters

		"""
		steadyWindsWanted = self.steadyWinds[self.steadyButtonsGroup.checkedId()]
		gustWindsWanted = self.gustWinds[self.gustButtonsGroup.checkedId()]
		# print(steadyWindsWanted, gustWindsWanted)
		self.aeroInstance.setWindModel(*steadyWindsWanted, gustWindsWanted)
		return
