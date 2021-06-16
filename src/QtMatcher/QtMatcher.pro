NAME         = QtMatcher
TARGET       = $${NAME}

QT           = core
QT          -= gui
QT          += network
QT          += sql
QT          += script
QT          += Essentials

load(qt_build_config)
load(qt_module)

INCLUDEPATH += $${PWD}/../../include/$${NAME}

HEADERS     += $${PWD}/../../include/$${NAME}/qtmatcher.h

SOURCES     += $${PWD}/nMatcher.cpp
SOURCES     += $${PWD}/nPatternExact.cpp
SOURCES     += $${PWD}/nPatternDice.cpp
SOURCES     += $${PWD}/nPatternCosine.cpp
SOURCES     += $${PWD}/nPatternJaccard.cpp
SOURCES     += $${PWD}/nPatternOverlap.cpp

OTHER_FILES += $${PWD}/../../include/$${NAME}/headers.pri

include ($${PWD}/../../doc/Qt/Qt.pri)

TRNAME       = $${NAME}
include ($${PWD}/../../Translations.pri)
