/*******************************************************************************
 * Copyright (C) 2017-2021 Theodore Chang
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 ******************************************************************************/
/**
 * @class Constraint
 * @brief A Constraint class.
 *
 * The Constraint class.
 *
 * @author tlc
 * @date 03/07/2017
 * @file Constraint.h
 * @addtogroup Constraint
 * @{
 */

#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include <Domain/Tag.h>
#include <Load/Amplitude/Amplitude.h>

class DomainBase;

class Constraint : public Tag {
protected:
	static const double multiplier;

	unsigned start_step = 0; /**< step tag */
	unsigned end_step = static_cast<unsigned>(-1);
	unsigned amplitude_tag = 0;

	uvec nodes; /**< node indices */
	uvec dofs;  /**< DoF indices */

	shared_ptr<Amplitude> magnitude;

	friend void set_constraint_multiplier(double);
public:
	const bool initialized = false;

	Constraint(unsigned, unsigned, unsigned, uvec&&, uvec&&);
	Constraint(const Constraint&) = delete;            // copy forbidden
	Constraint(Constraint&&) = delete;                 // move forbidden
	Constraint& operator=(const Constraint&) = delete; // assign forbidden
	Constraint& operator=(Constraint&&) = delete;      // assign forbidden

	virtual ~Constraint();

	virtual int initialize(const shared_ptr<DomainBase>&);

	virtual int process(const shared_ptr<DomainBase>&) = 0;

	void set_initialized(bool) const;

	void set_start_step(unsigned);
	[[nodiscard]] unsigned get_start_step() const;

	void set_end_step(unsigned);
	[[nodiscard]] unsigned get_end_step() const;

	[[nodiscard]] bool validate_step(const shared_ptr<DomainBase>&) const;

	// some constraint may manage state
	virtual void commit_status();
	virtual void clear_status();
	virtual void reset_status();
};

void set_constraint_multiplier(double);

#endif

//! @}
