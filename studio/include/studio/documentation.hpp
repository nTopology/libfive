/*
Studio: a simple GUI for the libfive CAD kernel
Copyright (C) 2017  Matt Keeter

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
*/
#pragma once

#include <QWidget>
#include <QLineEdit>
#include <QString>
#include <QMap>

namespace Studio {

typedef QMap<QString, QMap<QString, QString>> Documentation;

////////////////////////////////////////////////////////////////////////////////

class DocumentationPane : public QWidget
{
    Q_OBJECT
public:
    DocumentationPane(Documentation docs);

    /*  Shows the documentation pane and makes the search window focused */
    void show();


protected:
    /*
     *  Hides the completer popup in the close event; otherwise, if the
     *  DocumentationPane is closed by pressing the button in the corner,
     *  the completer popup remains awkwardly floating by itself.
     */
    void closeEvent(QCloseEvent* event) override;

    /*
     *  Special-case handling of escape (hides the window)
     */
    bool eventFilter(QObject* object, QEvent* event) override;

    QLineEdit* m_search;
};
} // namespace Studio
